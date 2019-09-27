################################################################################
# R utility functions accompanying:
# Tanentzap AJ, Igea J, Johnston MG, Larcombe MJ. 2019. 
# Does evolutionary history correlate with contemporary extinction risk by influencing range size dynamics?
# The American Naturalist
#
# No guarantees provided with this code.
# Prepared by AJ Tanentzap (27 Sep 2019, ajt65 // @ // cam.ac.uk)
#
# Code to run the analyses associated with the paper is in the file "R code to reproduce analyses.R"
# This file just contains utility functions, most of which are directly from https://github.com/Ax3man/phylopath
################################################################################


check_models_data_tree <- function(model_set, data, tree, na.rm) {
  var_names <- lapply(model_set, colnames)
  # Check whether all causal models have the same variables.
  if (length(model_set) > 1 &
      (stats::var(lengths(model_set)) != 0 |
       any(lengths(sapply(var_names[-1], setdiff, var_names[[1]])) != 0))) {
    stop('All causal models need to include the same variables. Your
       model set includes the following variables:\n',
         paste(sort(unique(unlist(var_names))), collapse = '\n'),
         call. = FALSE)
  }
  data <- data[, unique(unlist(var_names))]
  # We force all character columns to factors
  char_cols <- sapply(data, is.character)
  data[char_cols] <- lapply(data[char_cols], as.factor)
  # Check whether all factors have exactly two levels:
  f_cols <- which(sapply(data, is.factor))
  for (i in f_cols) {
    n_levels <- length(levels(data[[i]]))
    if (n_levels != 2) {
      stop("Variable '", names(data)[i], "' is expected to binary, but has ", n_levels, " levels.",
           .call = FALSE)
    }
  }
  # Check tree
  if (inherits(tree, 'multiPhylo')) {
    stop('You are passing several trees (in a `multiPhylo` object). Please only pass one `phylo` object.')
  }
  if (!inherits(tree, 'phylo')) {
    stop('The tree needs to be of class `phylo`.')
  }
  # Check NAs and if models and tree line up
  if (anyNA(data)) {
    if (na.rm) {
      NAs <- which(apply(data, 1, anyNA))
      message(length(NAs), ' rows were dropped because they contained NA values.')
      data <- data[-NAs, ]
    } else {
      stop('NA values were found in the variables of interest.', call. = FALSE)
    }
  }
  # Match the tree
  if (length(setdiff(rownames(data), tree$tip.label)) > 0) {
    stop('Make sure that species in your data have rownames that are exactly matched by name with tips in the tree.')
  }
  # Prune the tree
  if (length(tree$tip.label) > nrow(data)) {
    tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(data)))
    message('Pruned tree to drop species not included in dat.')
  }
  # Add names to the models, if they don't have them
  if (is.null(names(model_set))) {
    names(model_set) <- LETTERS[1:length(model_set)]
  }
  return(list(model_set = model_set, data = data, tree = tree))
}

find_consensus_order <- function(model_set) {
  # If the fully combined model is acyclic, then we use that.
  model_set_same_order <- lapply(model_set, function(x) {
    x[rownames(model_set[[1]]), colnames(model_set[[1]])]
  } )
  full_model <- sign(Reduce('+', model_set_same_order))
  if (ggm::isAcyclic(full_model)) {
    return(rownames(ggm::topSort(full_model)))
  }
  # Otherwise we find the most common orderings and use those.
  # Make sure all models are ordered:
  model_set <- lapply(model_set, ggm::topSort)
  vars <- lapply(model_set, row.names)
  combs <- as.data.frame(t(utils::combn(vars[[1]], 2)), stringsAsFactors = FALSE)
  names(combs) <- c('node1', 'node2')
  combs$count <- 0
  for (i in seq_along(vars)) {
    v <- apply(combs, 1, function(x) {
      which(vars[[i]] == x[1]) < which(vars[[i]] == x[2])
    } )
    combs$count <- combs$count + v
  }

  # If node1 is commonly ordered above node2, leave as is, otherwise swap them around
  combs <- dplyr::mutate_(combs,
                          tmp = ~node1,
                          node1 = ~ifelse(count > 0.5 * dplyr::n(), node1, node2),
                          node2 = ~ifelse(count > 0.5 * dplyr::n(), node2, tmp))
  # Now we order the nodes by how many nodes they are above, this should go from n:1
  combs <- dplyr::group_by_(combs, ~node1)
  combs <- dplyr::mutate_(combs, n = ~dplyr::n())
  combs <- dplyr::arrange_(combs, ~dplyr::desc(n))
  res <- unlist(c(unique(combs$node1), utils::tail(combs, 1)[, 'node2']))
  names(res) <- NULL
  res
}

set_to_formula <- function(x) {
  dep <- x[2]
  ind <- x[1]
  cond <- x[c(-1, -2)]

  stats::formula(paste(dep, paste(c(cond, ind), collapse = '+'), sep = '~'))
}

find_formulas <- function(d, order) {
  s <- ggm::basiSet(d)
  if (is.null(s)) {
    stop('One or some of your models are fully connected, and cannot be tested.')
  }
  s <- lapply(s, function(x) {
    # define whether there are existing paths between the two nodes in both directions.
    path1 <- !is.null(ggm::findPath(d, which(rownames(d) == x[1]), which(rownames(d) == x[2])))
    path2 <- !is.null(ggm::findPath(d, which(rownames(d) == x[2]), which(rownames(d) == x[1])))
    if (path1 & !path2) {
      # the first vertex is upstream, so we do not re-order
      return(x)
    }
    if ((path2 & !path1) | (path1 & path2)) {
      # these conditions should not occur, the first means basiSet is returning the wrong order,
      # the second should only occur if there are cycles.
      stop('If you get this error, please contact the maintainer.')
    }
    if (!path1 & !path2) {
      # check whether the order is according to `order`
      if (which(order == x[1]) < which(order == x[2])) {
        return(x)
      } else {
        return(c(x[2], x[1], x[-(1:2)]))
      }
    }
  } )
  lapply(s, set_to_formula)
}

C_stat <- function(ps) -2 * sum(log(ps))

C_p <- function(C, k) 1 - stats::pchisq(C, 2 * k)

CICc <- function(C, q, n) C + 2 * q * (n / (n - 1 - q))

l <- function(dCICc) exp(-0.5 * dCICc)

w <- function(l) l / sum(l)

phylo_g_lm2 <- function(formula, data, tree, model, method, boot = 0, ...) {
  # we capture the dots, because we need to match the names to either phylolm or phylolm
  dots <- list(...)
  dots_glm <- dots[names(dots) %in% names(formals(phylolm::phyloglm))]
  dots_lm <- dots[names(dots) %in% names(formals(phylolm::phylolm))]
  if (length(intersect(names(dots_glm), names(dots_lm))) != length(dots)) {
    warning("Some arguments in ... are not recognized.", call. = FALSE)
  }
  # we capture the first argument in the formula, to check whether it is binary
  x_var <- data[[all.vars(formula)[1]]]
  if (is.factor(x_var)) {
    # phyloglm need binary variables as 0,1 but I use factors
    data[all.vars(formula)[1]] <- as.numeric(x_var) - 1
    fun <- phylolm::phyloglm
    args <- c(list(formula = formula, data = data, phy = tree, method = method, boot = boot, log.alpha.bound = 20, btol = 20),
              dots_glm)
  } else {
    fun <- phylolm::phylolm
    args <- c(list(formula = formula, data = data, phy = tree, model = model, boot = boot),
              dots_glm)
  }
  res <- do.call(quiet_safely(fun), args)
  # Remove the call, since quiet_safely messes it up and it's annoying in printing
  res$result$call <- NULL

  return(res)
}

get_p <- function(m) {
  s <- stats::coef(summary(m))
  return(s[nrow(s), 'p.value'])
}

get_est <- function(m) {
  stats::coef(m)[-1]
}

get_se <- function(m) {
  stats::coef(summary(m))[-1, 'StdErr']
}

get_lower <- function(m) {
  s <- stats::coef(summary(m))
  if ('lowerbootCI' %in% colnames(s)) {
    r <- s[-1, 'lowerbootCI']
  } else {
    r <- NA
  }
  return(r)
}

get_upper <- function(m) {
  s <- stats::coef(summary(m))
  if ('upperbootCI' %in% colnames(s)) {
    r <- s[-1, 'upperbootCI']
  } else {
    r <- NA
  }
  return(r)
}

get_phylo_param <- function(m) {
  r <- m$optpar
  if (is.null(r)) r <- NA
  return(r)
}

adjust_layout <- function(l, rotation, flip_x, flip_y) {
  rotation <- rotation * (2 * pi / 360)
  R <- matrix(c(cos(rotation), sin(rotation), -sin(rotation), cos(rotation)), nrow = 2)
  l[c('x', 'y')] <- as.matrix(l[c('x', 'y')]) %*% R
  if (flip_x) {
    l$x <- -l$x
  }
  if (flip_y) {
    l$y <- -l$y
  }
  return(l)
}

combine_with_labels <- function(l, labels) {
  if (is.null(labels)) {
    return(l)
  }
  if (is.null(names(labels))) {
    stop('labels must be a named vector.', call. = FALSE)
  }
  if (length(setdiff(l$name, names(labels))) > 0) {
    stop('Some nodes are missing from labels.', call. = FALSE)
  }
  l$name <- factor(l$name, names(labels), labels)
  class(l) <- c("layout_igraph", "layout_ggraph", "data.frame")
  return(l)
}

quiet_safely <- function(.f) {
  capture_all <- function(expr)
  {
    warn_vec <- NULL
    w.handler <- function(w){ # warning handler
      warn_vec <<- c(warn_vec, w$message)
      invokeRestart("muffleWarning")
    }
    r <- list(result = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
              warning = warn_vec)
    if (inherits(r$result, 'error')) {
      return(list(result = NULL, error = r$result$message, warning = r$warning))
    } else {
    return(list(result = r$result, error = NULL, warning = r$warning))
    }
  }
  function(...) capture_all(.f(...))
}

combine_dots <- function(old_dots, ...) {
  new_dots <- list(...)
  c(new_dots, old_dots[!(names(old_dots) %in% names(new_dots))])
}



#' Compare causal models in a phylogenetic context.
#'
#' Continuous variables are modeled using [phylolm::phylolm], while binary
#' traits are modeled using [phylolm::phyloglm].
#'
#' @param model_set A list of directed acyclic graphs. These are matrices,
#'   typically created with \code{define_model_set}.
#' @param data A \code{data.frame} with data. If you have binary variables, make
#'   sure they are either character values or factors!
#' @param tree A phylogenetic tree of class \code{phylo}.
#' @param model The evolutionary model used for the regressions on continuous
#'   variables. See [phylolm::phylolm] for options and details. Defaults to
#'   Pagel's lambda model
#' @param method The estimation method for the binary models. See
#'   [phylolm::phylolm] for options and details. Defaults to logistic MPLE.
#' @param order Causal order of the included variable, given as a character
#'   vector. This is used to determine which variable should be the dependent
#'   in the dsep regression equations. If left unspecified, the order will be
#'   automatically determined. If the combination of all included models is
#'   itself a DAG, then the ordering of that full model is used. Otherwise,
#'   the most common ordering between each pair of variables is used to create
#'   a general ordering.
#' @param parallel An optional vector containing the virtual connection
#'   process type for running the chains in parallel (such as \code{"SOCK"}).
#'   A cluster is create using the \code{parallel} package.
#' @param na.rm Should rows that contain missing values be dropped from the data
#'   as necessary (with a message)?
#' @param ...
#'   Arguments passed on to `phylolm`:
#'
#'   `lower.bound`: optional lower bound for the optimization of the phylogenetic model parameter.
#'
#'   `upper.bound`: optional upper bound for the optimization of the phylogenetic model parameter.
#'
#'   `starting.value`: optional starting value for the optimization of the phylogenetic model parameter.
#'
#'   `measurement_error`: a logical value indicating whether there is measurement error sigma2_error (see Details).
#'
#'   Arguments passed on to `phyloglm`:
#'
#'   `btol`: bound on the linear predictor to bound the searching space.
#'
#'   `log.alpha.bound`: bound for the log of the parameter alpha.
#'
#'   `start.beta`: starting values for beta coefficients.
#'
#'   `start.alpha`: starting values for alpha (phylogenetic correlation).
#'
#'
#' @return A phylopath object, with the following components:
#'  \describe{
#'   \item{d_sep}{for each model a table with separation statements and statistics.}
#'   \item{model_set}{the DAGs}
#'   \item{data}{the supplied data}
#'   \item{tree}{the supplied tree}
#'   \item{model}{the employed model of evolution in `phylolm`}
#'   \item{method}{the employed method in `phyloglm`}
#'   \item{dots}{any additional arguments given, these are passed on to downstream functions}
#'   \item{warnings}{any warnings generated by the models}
#'   }
#' @export
#' @examples
#'   #see vignette('intro_to_phylopath') for more details
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'
#'   # Printing p gives some general information:
#'   p
#'   # And the summary gives statistics to compare the models:
#'   summary(p)
#'
phylo_path2 <- function(model_set, data, tree, model = 'lambda', method = 'logistic_MPLE',
                       order = NULL, parallel = NULL, na.rm = TRUE, ...) {
  # Always coerce to data.frame, as tibbles and data.tables do NOT play nice.
  data <- as.data.frame(data)
  tmp <- check_models_data_tree(model_set, data, tree, na.rm)
  model_set <- tmp$model_set
  data <- tmp$data
  tree <- tmp$tree

  if (is.null(order)) {
    order <- find_consensus_order(model_set)
  }
  formulas <- lapply(model_set, find_formulas, order)
  formulas <- purrr::map(formulas,
                         ~purrr::map(.x, ~{attr(., ".Environment") <- NULL; .}))
  f_list <- unique(unlist(formulas))
  if (!is.null(parallel)) {
    cl <- parallel::makeCluster(min(c(parallel::detectCores() - 1,
                                      length(f_list))),
                                parallel)
    parallel::clusterExport(cl, list('phylo_g_lm2'), environment())
    on.exit(parallel::stopCluster(cl))
  } else {
    cl <- NULL
  }
  dsep_models_runs <- pbapply::pblapply(
    f_list,
    function(x, data, tree, model, method, ...) {
      phylo_g_lm2(x, data, tree, model, method, ...)
    },
    data = data, tree = tree, model = model, method = method, cl = cl)
  # Produce appropriate error if needed
  errors <- purrr::map(dsep_models_runs, 'error')
  purrr::map2(errors, f_list,
              ~if(!is.null(.x))
                stop(paste('Fitting the following model:\n   ',
                           Reduce(paste, deparse(.y)),
                           '\nproduced this error:\n   ', .x),
                     call. = FALSE))
  # Collect warnings as well, but save those for later.
  warnings <- purrr::map(dsep_models_runs, 'warning')
  warnings <- purrr::map2(warnings, f_list,
                          ~if(!is.null(.x))
                             paste('Fitting the following model:\n   ',
                                       Reduce(paste, deparse(.y)),
                                       '\nproduced this/these warning(s):\n   ', .x))
  warnings <- warnings(!sapply(warnings, is.null))
  if (length(warnings) > 1) {
    warning('Some models produced warnings. Use `show_warnings()` to view them.')
  }

  # Collect models.
  dsep_models <- purrr::map(dsep_models_runs, 'result')
  dsep_models <- purrr::map(formulas, ~dsep_models[match(.x, f_list)])

  d_sep <- purrr::map2(
    formulas,
    dsep_models,
    ~dplyr::data_frame(
      d_sep = as.character(.x),
      p = purrr::map_dbl(.y, get_p),
      phylo_par = purrr::map_dbl(.y, get_phylo_param),
      model = .y
    )
  )

  out <- list(d_sep = d_sep, model_set = model_set, data = data, tree = tree,
              model = model, method = method, dots = list(...), warnings = warnings)
  class(out) <- 'phylopath'
  return(out)
}

#' @export
summary.phylopath <- function(object, ...) {
  phylopath <- object
  stopifnot(inherits(phylopath, 'phylopath'))
  k <- sapply(phylopath$d_sep, nrow)
  q <- sapply(phylopath$model_set, function(m) nrow(m) + sum(m))
  C <- sapply(phylopath$d_sep, function(x) C_stat(x$p))
  p <- C_p(C, k)
  IC <- CICc(C, q, nrow(phylopath$data))

  d <- data.frame(model = names(phylopath$model_set), k = k, q = q, C = C, p = p,
                  CICc = IC, stringsAsFactors = FALSE)
  d <- d[order(d$CICc), ]
  d$delta_CICc <- d$CICc - d$CICc[1]
  d$l <- l(d$delta_CICc)
  d$w <- w(d$l)
  class(d) <- c('phylopath_summary', 'data.frame')
  return(d)
}

#' Extract and estimate the best supported model from a phylogenetic path
#' analysis.
#'
#' @param phylopath An object of class \code{phylopath}.
#' @param ... Arguments to pass to [phylolm::phylolm] and [phylolm::phyloglm]. If you specified
#'   options in the original [phylo_path] call you don't need to specify them again.
#'
#' @return An object of class \code{fitted_DAG}.
#' @export
#'
#' @examples
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   best_model <- best(p)
#'   # Print the best model to see coefficients, se and ci:
#'   best_model
#'   # Plot to show the weighted graph:
#'   plot(best_model)
#'
best <- function(phylopath, ...) {
  stopifnot(inherits(phylopath, 'phylopath'))
  dots <- combine_dots(phylopath$dots, ...)

  b <- summary(phylopath)[1, 'model']
  best_model <- phylopath$model_set[[b]]
  do.call(
    est_DAG,
    c(list(best_model, phylopath$data, phylopath$tree, phylopath$model, phylopath$method), dots)
  )
}

#' Extract and estimate an arbitrary model from a phylogenetic path analysis.
#'
#' @param phylopath An object of class \code{phylopath}.
#' @param choice A character string of the name of the model to be chosen, or
#'   the index in \code{model_set}.
#' @param ... Arguments to pass to [phylolm::phylolm] and [phylolm::phyloglm]. If you specified
#'   options in the original [phylo_path] call you don't need to specify them again.
#'
#' @return An object of class \code{fitted_DAG}.
#' @export
#'
#' @examples
#'   candidates <- list(A = DAG(LS ~ BM, NL ~ BM, DD ~ NL),
#'                      B = DAG(LS ~ BM, NL ~ LS, DD ~ NL))
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   my_model <- choice(p, "B")
#'   # Print the best model to see coefficients, se and ci:
#'   my_model
#'   # Plot to show the weighted graph:
#'   plot(my_model)
#'
choice <- function(phylopath, choice, ...) {
  stopifnot(inherits(phylopath, 'phylopath'))
  dots <- combine_dots(phylopath$dots, ...)

  do.call(
    est_DAG,
    c(list(phylopath$model_set[[choice]], phylopath$data, phylopath$tree, phylopath$model,
           phylopath$method), dots)
  )
}

#' Extract and average the best supported models from a phylogenetic path
#' analysis.
#'
#' @param phylopath An object of class `phylopath`.
#' @param cut_off The CICc cut-off used to select the best models. Use
#'   `Inf` to average over all models. Use the [best()] function to
#'   only use the top model, or [choice()] to select any single model.
#' @param ... Arguments to pass to [phylolm::phylolm] and [phylolm::phyloglm]. If you specified
#'   options in the original [phylo_path] call you don't need to specify them again.
#'
#' @inheritParams average_DAGs
#'
#' @return An object of class `fitted_DAG`.
#' @export
#'
#' @examples
#'   candidates <- list(
#'     A = DAG(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ BM + NL),
#'     B = DAG(LS ~ BM, NL ~ BM + RS, DD ~ NL)
#'   )
#'   p <- phylo_path(candidates, rhino, rhino_tree)
#'   summary(p)
#'
#'   # Models A and B have similar support, so we may decide to take
#'   # their average.
#'
#'   avg_model <- average(p)
#'   # Print the average model to see coefficients, se and ci:
#'   avg_model
#'
#'   \dontrun{
#'   # Plot to show the weighted graph:
#'   plot(avg_model)
#'
#'   # One can see that an averaged model is not necessarily a DAG itself.
#'   # This model actually has a path in two directions.
#'
#'   # Note that coefficients that only occur in one of the models become much
#'   # smaller when we use full averaging:
#'
#'   coef_plot(avg_model)
#'   coef_plot(average(p, method = 'full'))
#'   }
#'
average <- function(phylopath, cut_off = 2, avg_method = 'conditional', ...) {
  stopifnot(inherits(phylopath, 'phylopath'))
  dots <- combine_dots(phylopath$dots, ...)

  d <- summary(phylopath)
  b <- d[d$delta_CICc < cut_off, ]

  best_models <- lapply(
    phylopath$model_set[b$model],
    function(x) {
      do.call(
        est_DAG,
        c(list(DAG = x, data = phylopath$data, tree = phylopath$tree, model = phylopath$model,
               method = phylopath$method), dots)
      )
    }
  )
  average <- average_DAGs(best_models, b$w, avg_method)
  class(average$coef) <- c('matrix', 'DAG')
  return(average)
}


phylo_path2 <- function(model_set, data, tree, model = 'lambda', method = 'logistic_MPLE',
                        order = NULL, parallel = NULL, na.rm = TRUE, ...) {
  # Always coerce to data.frame, as tibbles and data.tables do NOT play nice.
  data <- as.data.frame(data)
  tmp <- check_models_data_tree(model_set, data, tree, na.rm)
  model_set <- tmp$model_set
  data <- tmp$data
  tree <- tmp$tree

  if (is.null(order)) {
    order <- find_consensus_order(model_set)
  }
  formulas <- lapply(model_set, find_formulas, order)
  formulas <- purrr::map(formulas,
                         ~purrr::map(.x, ~{attr(., ".Environment") <- NULL; .}))
  f_list <- unique(unlist(formulas))
  if (!is.null(parallel)) {
    cl <- parallel::makeCluster(min(c(parallel::detectCores() - 1,
                                      length(f_list))),
                                parallel)
    parallel::clusterExport(cl, list('phylo_g_lm2'), environment())
    on.exit(parallel::stopCluster(cl))
  } else {
    cl <- NULL
  }
  dsep_models_runs <- pbapply::pblapply(
    f_list,
    function(x, data, tree, model, method, ...) {
      phylo_g_lm2(x, data, tree, model, method, ...)
    },
    data = data, tree = tree, model = model, method = method, cl = cl)
  # Produce appropriate error if needed
  errors <- purrr::map(dsep_models_runs, 'error')
  purrr::map2(errors, f_list,
              ~if(!is.null(.x))
                stop(paste('Fitting the following model:\n   ',
                           Reduce(paste, deparse(.y)),
                           '\nproduced this error:\n   ', .x),
                     call. = FALSE))
  # Collect warnings as well, but save those for later.
  warnings <- purrr::map(dsep_models_runs, 'warning')
  warnings <- purrr::map2(warnings, f_list,
                          ~if(!is.null(.x))
                             paste('Fitting the following model:\n   ',
                                       Reduce(paste, deparse(.y)),
                                       '\nproduced this/these warning(s):\n   ', .x))
  warnings <- warnings(!sapply(warnings, is.null))
  if (length(warnings) > 1) {
    warning('Some models produced warnings. Use `show_warnings()` to view them.')
  }

  # Collect models.
  dsep_models <- purrr::map(dsep_models_runs, 'result')
  dsep_models <- purrr::map(formulas, ~dsep_models[match(.x, f_list)])

  d_sep <- purrr::map2(
    formulas,
    dsep_models,
    ~dplyr::data_frame(
      d_sep = as.character(.x),
      p = purrr::map_dbl(.y, get_p),
      phylo_par = purrr::map_dbl(.y, get_phylo_param),
      model = .y
    )
  )

  out <- list(d_sep = d_sep, model_set = model_set, data = data, tree = tree,
              model = model, method = method, dots = list(...), warnings = warnings)
  class(out) <- 'phylopath'
  return(out)
}
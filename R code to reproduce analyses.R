################################################################################
# R code accompanying:
# Tanentzap AJ, Igea J, Johnston MG, Larcombe MJ. 2019. 
# Does evolutionary history correlate with contemporary extinction risk by influencing range size dynamics?
# The American Naturalist
#
# No guarantees provided with this code.
# Prepared by AJ Tanentzap (18 Sep 2019, ajt65 // @ // cam.ac.uk)
#
# Code includes the following sections:
# 1) load packages and relevant functions
# 2) load data for genus-level analysis
# 3) test whether threat status associated with diversification and age for genus-level analysis
# 4) are we taking a bias part of the tol for genus-level analysis?
# 5) how similar are genus- and species-level ages?
# 6) load data for species-level analysis
# 7) test association among age, threat status, and range in conifers
# 8) sister species analysis in conifers
# 9) rerun analysis at genus-level in conifers
# 10) check bias in the conifer dataset
# 11) test association among age, threat status, and range in palms
# 12) sister species analysis in palms
# 13) rerun analysis at genus-level in palms
# 14) check bias in the palm dataset
# 15) power analysis for the genus-level results
# 16) sensitivity analysis for sampling coverage

# Data accompanying code available in "./data/" folder of the GitHub repo https://github.com/atanzap/age-extinction-plants
# Code assumes working directory set to that containing data files.

################################################################################
###### 1) LOAD PACKAGES AND RELEVANT FUNCTIONS
################################################################################
library(ape)
library(Taxonstand)
library(geiger)
library(taxonlookup)
library(phytools)
library(caper)
library(rredlist)
library(phyndr)
library(MonoPhy)
library(picante)
library(phylolm)
library(car)
library(brms)
library(MCMCglmm)
library(phylopath)
library(nlme)
library(rr2)

#function to calculate species ages given that tree is a ape phylo object
calculate_tip_ages<-function(tree){
  node.age(tree)->phy.age
  cbind(phy.age$edge,phy.age$age, tree$edge.length)->BL.position
  max(phy.age$age)-BL.position[,3]->dist.tip
  cbind(BL.position,dist.tip)->BL.positions
  BL.positions[,5]+BL.positions[,4]->ages
  cbind(BL.positions,ages)->BL.positions
  as.data.frame(BL.positions)->node.ages
  names(node.ages)<-c("parental.node","daughter.node","dist.root","BL","dist.tip","mrca.age")
  ## node.ages is a data frame listing as variables the identity of parental and
  #daughter nodes, the distance from the root and from the present of each node,
  #the branch length and the age of the most recent common ancestor
  node.ages[node.ages[,2]<length(tree$tip)+1,]->species.ages
  row.names(species.ages)<-tree$tip[species.ages$daughter.node]
  ## species ages is node.ages data frame reduced to the tips (species)
  species.ages<-species.ages[order(row.names(species.ages)),]
  output.table<-as.data.frame(cbind(row.names(species.ages),species.ages$mrca.age))
  colnames(output.table)<-c('tip','tip.age')
  return(output.table)
}

# function to extract characteristics for data frame of sisters
# x = sister species to analyse, dataf = dataframe with their characteristics, sampling = whether to permutate sisters (i.e. randomly sample difference)
sis_chars <- function(x, dataf, sampling){ 
                          min_diverg <- min(dataf[which(rownames(dataf) %in% as.character(unlist(x[1:2]))),]$tip.age)
						              ifelse(sampling==T,ord<-sample(1:2,2),ord<-1:2)
                          diff1 <- diff(log(dataf[match(c(as.character(unlist(x[1:2]))[which(x[3:4] == 1)], c(as.character(unlist(x[1:2]))[which(x[3:4] == 0)])),rownames(dataf)),]$global.count)[ord])
                          return(c(x[1:2],min_diverg,diff1))
                        }

# functions to ensure phyloglm seeing changes to log.alpha.bound and btol when used in phylo_path
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





################################################################################
###### 2) LOAD DATA FOR GENUS-LEVEL ANALYSIS
################################################################################
# read IUCN data
categories.species<-read.csv('export-93437.csv',header=T)
categories.species$genus<-sapply(as.character(categories.species$Genus),function(x) strsplit(x,split=' ')[[1]][[1]])
categories.species$species<-sapply(as.character(categories.species$Species),function(x) strsplit(x,split=' ')[[1]][[1]])
# run through TPL
IUCN_TPL<-TPL(genus=categories.species$genus,species=categories.species$species)
IUCN_TPL_risk<-cbind(IUCN_TPL,categories.species[match(with(IUCN_TPL,paste(Genus,Species,sep=" ")),with(categories.species,paste(genus,species, sep= " "))),])
IUCN_TPL_risk[,'Merged']<-paste(IUCN_TPL_risk$New.Genus,IUCN_TPL_risk$New.Species, sep = '_')
IUCN_TPL_risk <- IUCN_TPL_risk[-which(IUCN_TPL_risk$Taxonomic.status == ""),]
IUCN_dups <- table(IUCN_TPL_risk$Merged)
IUCN_TPL_riskunique <- IUCN_TPL_risk[-which(IUCN_TPL_risk$Merged %in% names(IUCN_dups[which(IUCN_dups > 1)]) & (with(IUCN_TPL_risk,paste(Genus,Species,sep="_")) %in% IUCN_TPL_risk$Merged == F)),]
# run taxonloookup
IUCN_TPL_Lookup <- lookup_table(as.character(IUCN_TPL_riskunique$Merged),by_species = TRUE,include_counts = TRUE, family.tax="tpl")
rownames(IUCN_TPL_Lookup) <- NULL
IUCN_TPL_Lookup <- IUCN_TPL_Lookup[!duplicated(IUCN_TPL_Lookup$genus),]

# process Qian tree
PhytoPhylo <- read.tree('PhytoPhylo.nex')
PhytoPhylo_names <- sapply(PhytoPhylo$tip.label,strsplit,split='_')

# tidy hybrid names
PhytoPhylo_names[which(sapply(PhytoPhylo_names, length) == 3)] <- lapply(PhytoPhylo_names[which(sapply(PhytoPhylo_names, length) == 3)], function(x){c(x[1],x[3])})

# run tree through TPL
PhytoPhylo_TPL <- TPL(genus=as.vector(sapply(PhytoPhylo_names, function(x){x[1]})),species=as.vector(sapply(PhytoPhylo_names, function(x){x[2]})))

# tidy variety / subspecies names
PhytoPhylo_names_infras <- PhytoPhylo_names[which(sapply(PhytoPhylo_names, length) == 4)]
PhytoPhylo_TPL[which(sapply(PhytoPhylo_names, length) == 4),] <- TPL(genus=as.vector(sapply(PhytoPhylo_names_infras, function(x){x[1]})),species=as.vector(sapply(PhytoPhylo_names_infras, function(x){x[2]})), infrasp=as.vector(sapply(PhytoPhylo_names_infras, function(x){x[4]})))
         
# rename phylo object
PhytoPhylo$tip.label <- with(PhytoPhylo_TPL, paste(New.Genus,New.Species,sep='_'))

#checking for monophyly
monoreport <- AssessMonophyly(PhytoPhylo)
outliertips <- GetOutlierTips(monoreport)
outliertips <- unlist(outliertips)
outliertips <- unname(outliertips)
#drop the outliers from tree
PhytoPhylo_nooutliers<-drop.tip(PhytoPhylo,outliertips)

#genus level tree, get genera only
generaPhyndr <- unique(sapply(strsplit(as.character(PhytoPhylo_nooutliers$tip.label),'_'),function(x) x[1]))
generaPhyndr.2 <- sapply(generaPhyndr, function(x) x<- paste(x,'_',sep=''))
ii <- sapply(generaPhyndr.2,function(x,y) grep(x,y)[1],y=PhytoPhylo_nooutliers$tip.label)
#drop all but one of every genera
PhytoPhylo_GenusTree <- drop.tip(PhytoPhylo_nooutliers,setdiff(PhytoPhylo_nooutliers$tip.label,PhytoPhylo_nooutliers$tip.label[ii]))
PhytoPhylo_GenusTree <- drop.tip(PhytoPhylo_GenusTree,which(duplicated(PhytoPhylo_GenusTree$tip.label)))
#remove full species name and leave only genera as tip label
PhytoPhylo_GenusTree$tip.label <- sapply(strsplit(PhytoPhylo_GenusTree$tip.label,"_"),function(x) x[1])


# calculate tip ages
ages1 <- calculate_tip_ages(PhytoPhylo_GenusTree)
ages1[,2] <- as.numeric(as.character(ages1[,2]))


# for each genus calculate diversification
IUCN_ages <- merge(ages1,IUCN_TPL_Lookup,by.x="tip",by.y="genus",all.y=T)
IUCN_ages$r_hat0 <- bd.ms(time=IUCN_ages$tip.age, n=IUCN_ages$number.of.accepted.species, crown=F, epsilon = 0)
IUCN_ages$r_hat50 <- bd.ms(time=IUCN_ages$tip.age, n=IUCN_ages$number.of.accepted.species, crown=F, epsilon = 0.5)
IUCN_ages$r_hat90 <- bd.ms(time=IUCN_ages$tip.age, n=IUCN_ages$number.of.accepted.species, crown=F, epsilon = 0.9)


# find family-level sampling coverage in genus tree
plant_lookup_table <- plant_lookup(include_counts = TRUE, family.tax="tpl")
genera_in_families <- table(plant_lookup_table[plant_lookup_table$number.of.accepted.species >0,]$family)
IUCN_ages$number.of.genera <- as.numeric(genera_in_families[match(IUCN_ages$family, names(genera_in_families))])
plant_lookup_table_phylo <- plant_lookup_table[match(PhytoPhylo_GenusTree$tip.label,plant_lookup_table$genus),]
genera_per_family_phylo <- table(plant_lookup_table_phylo[which(plant_lookup_table_phylo$number.of.accepted.species>0),'family'])
IUCN_ages$prop.genera.sampled <- as.numeric(genera_per_family_phylo[match(IUCN_ages$family, names(genera_per_family_phylo))])/IUCN_ages$number.of.genera
IUCN_ages$prop.genera.sampled[is.na(IUCN_ages$prop.genera.sampled)] <- 0

# calculate proportion of assessments that are threatened, which accounts for biases in surveying threatened taxa
IUCN_TPL_riskunique_drop <- IUCN_TPL_riskunique[which( is.na(IUCN_TPL_riskunique$Red.List.status) == F & IUCN_TPL_riskunique$Red.List.status %in% c('DD','EW','EX') ==F ),]
threats_per_genus <- with(IUCN_TPL_riskunique_drop, tapply(Red.List.status %in% c("CR","EN","VU"), New.Genus, sum))
threats_per_genus_criterionB <- tapply(as.numeric(IUCN_TPL_riskunique_drop$Red.List.status %in% c("CR","EN","VU") & ( grepl("B", IUCN_TPL_riskunique_drop$Red.List.criteria) |
                                                                                                                      grepl("D2", IUCN_TPL_riskunique_drop$Red.List.criteria)   |
                                                                                                                      IUCN_TPL_riskunique_drop$Red.List.criteria == 'D1+2') ),
                                IUCN_TPL_riskunique_drop$New.Genus, sum)
assessments_per_genus <- table(IUCN_TPL_riskunique_drop$New.Genus)
IUCN_ages$n_threat <- as.numeric(threats_per_genus[match(IUCN_ages$tip,names(threats_per_genus))])
IUCN_ages$n_threat_criterionB <- as.numeric(threats_per_genus_criterionB[match(IUCN_ages$tip,names(threats_per_genus_criterionB))])
IUCN_ages$n_assess <- as.numeric(assessments_per_genus[match(IUCN_ages$tip,names(assessments_per_genus))])


# drop genera from families with poor sampling coverage, are very young, monotypic (so no r_hat), or lack any assessments
IUCN_ages_sub <- IUCN_ages[which(IUCN_ages$prop.genera.sampled >= .6 & IUCN_ages$r_hat0 != 0 & IUCN_ages$n_assess/IUCN_ages$number.of.accepted.species >= 0.2 & IUCN_ages$number.of.accepted.species > 1),]
PhytoPhylo_IUCNsub <- drop.tip(PhytoPhylo_GenusTree,setdiff(PhytoPhylo_GenusTree$tip.label,IUCN_ages_sub$tip))
PhytoPhylo_IUCNsub$node.label <- NULL
IUCN_ages_sub <- IUCN_ages_sub[match(PhytoPhylo_IUCNsub$tip.label, IUCN_ages_sub$tip),]
#drop ferns for now because there are too few
IUCN_ages_sub <- IUCN_ages_sub[IUCN_ages_sub$group != 'Pteridophytes',]





################################################################################
###### 3) TEST WHETHER THREAT STATUS ASSOCIATED WITH DIVERSIFICATION AND AGE FOR GENUS-LEVEL ANALYSIS
################################################################################
IUCN_PhytoPhylo_merged <- comparative.data(phy=PhytoPhylo_IUCNsub,data=IUCN_ages_sub,names.col=tip,vcv=T)
gls_wts <- sqrt(IUCN_ages_sub$n_assess/IUCN_ages_sub$number.of.accepted.species)
# run PGLS weighted by the number of species that are assessed
mod1b <- gls(I(n_threat_criterionB/n_assess) ~ log(r_hat0), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod2b <- gls(I(n_threat_criterionB/n_assess) ~ log(r_hat50), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod3b <- gls(I(n_threat_criterionB/n_assess) ~ log(r_hat90), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
#mod1b_asin <- gls(asin(I(n_threat_criterionB/n_assess)) ~ log(r_hat0), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
#mod2b_asin <- gls(asin(I(n_threat_criterionB/n_assess)) ~ log(r_hat50), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
#mod3b_asin <- gls(asin(I(n_threat_criterionB/n_assess)) ~ log(r_hat90), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
R2.pred(mod2b)

# repeat age analysis but don't drop monotypic genera
IUCN_ages_sub2 <- IUCN_ages[which(IUCN_ages$prop.genera.sampled >= .6 & IUCN_ages$n_assess/IUCN_ages$number.of.accepted.species >= 0.2 & IUCN_ages$group != 'Pteridophytes' & IUCN_ages$number.of.accepted.species > 0),]
PhytoPhylo_IUCNsub2 <- drop.tip(PhytoPhylo_GenusTree,setdiff(PhytoPhylo_GenusTree$tip.label,IUCN_ages_sub2$tip))
PhytoPhylo_IUCNsub2$node.label <- NULL
IUCN_ages_sub2 <- IUCN_ages_sub2[match(PhytoPhylo_IUCNsub2$tip.label, IUCN_ages_sub2$tip),]
IUCN_PhytoPhylo_merged2 <- comparative.data(phy=PhytoPhylo_IUCNsub2,data=IUCN_ages_sub2,names.col=tip,vcv=T)
gls_wts2 <- sqrt(IUCN_ages_sub2$n_assess/IUCN_ages_sub2$number.of.accepted.species)
mod4b <- gls(I(n_threat_criterionB/n_assess) ~ log(tip.age), data=IUCN_ages_sub2, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged2$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts2))
#mod4b_asin <- gls(asin(I(n_threat_criterionB/n_assess)) ~ log(tip.age), data=IUCN_ages_sub2, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged2$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts2))
R2.pred(mod4b)





################################################################################
###### 4) ARE WE TAKING A BIAS PART OF THE TOL FOR GENUS-LEVEL ANALYSIS?
################################################################################
# age distribution across tree
all_ages <- merge(ages1,plant_lookup_table,by.x="tip",by.y="genus",all.x=T)
all_ages$number.of.genera <- as.numeric(genera_in_families[match(all_ages$family, names(genera_in_families))])
all_ages$prop.genera.sampled <- as.numeric(genera_per_family_phylo[match(all_ages$family, names(genera_per_family_phylo))])/all_ages$number.of.genera
all_ages$r_hat0 <- bd.ms(time=all_ages$tip.age, n=all_ages$number.of.accepted.species, crown=F, epsilon = 0)
leveneTest(c(log(IUCN_ages_sub$tip.age),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$number.of.accepted.species > 1 & all_ages$r_hat0 != 0),]$tip.age)),
           c(rep("data",dim(IUCN_ages_sub)[1]),rep("all",dim(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$number.of.accepted.species > 1 & all_ages$r_hat0 != 0),])[1])))
t.test(log(IUCN_ages_sub$tip.age),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$number.of.accepted.species > 1 & all_ages$r_hat0 != 0),]$tip.age), var.equal=T)
mean(log(IUCN_ages_sub$tip.age))
sd(log(IUCN_ages_sub$tip.age))/sqrt(length(log(IUCN_ages_sub$tip.age)))
mean(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$number.of.accepted.species > 1 & all_ages$r_hat0 != 0),]$tip.age))
sd(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$number.of.accepted.species > 1 & all_ages$r_hat0 != 0),]$tip.age))/sqrt(length(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$number.of.accepted.species > 1 & all_ages$r_hat0 != 0),]$tip.age)))

# diversity distribution across tree
leveneTest(c(log(IUCN_ages_sub$number.of.accepted.species),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species)),
           c(rep("data",dim(IUCN_ages_sub)[1]),rep("all",dim(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$tip.age != 0 & all_ages$number.of.accepted.species > 1),])[1])))
wilcox.test(log(IUCN_ages_sub$number.of.accepted.species),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species), var.equal=F)
mean(log(IUCN_ages_sub$number.of.accepted.species))
sd(log(IUCN_ages_sub$number.of.accepted.species))/sqrt(length(log(IUCN_ages_sub$number.of.accepted.species)))
mean(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species))
sd(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species))/sqrt(length(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species)))

# rate distribution across tree
all_ages$r_hat50 <- bd.ms(time=all_ages$tip.age, n=all_ages$number.of.accepted.species, crown=F, epsilon = 0.50)
leveneTest(c(log(IUCN_ages_sub$r_hat50),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat50 != 0 & all_ages$number.of.accepted.species > 1),]$r_hat50)),
           c(rep("data",dim(IUCN_ages_sub)[1]),rep("all",dim(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat50 != 0 & all_ages$number.of.accepted.species > 1),])[1])))
t.test(log(IUCN_ages_sub$r_hat50),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat50 != 0 & all_ages$number.of.accepted.species > 1),]$r_hat50), var.equal=T)
mean(log(IUCN_ages_sub$r_hat50))
sd(log(IUCN_ages_sub$r_hat50))/sqrt(length(log(IUCN_ages_sub$r_hat50)))
mean(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat50 != 0 & all_ages$number.of.accepted.species > 1),]$r_hat50))
sd(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat50 != 0 & all_ages$number.of.accepted.species > 1),]$r_hat50))/sqrt(length(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat50 != 0 & all_ages$number.of.accepted.species > 1),]$r_hat50)))

# percent of threatened species
IUCN_allgens <- as.data.frame(cbind(threats_per_genus_criterionB, assessments_per_genus))
IUCN_allgens$genus <- as.character(rownames(IUCN_allgens))
rownames(IUCN_allgens) <- NULL
IUCN_allgens <- merge(IUCN_allgens,plant_lookup_table,all.x=T,by='genus')
IUCN_allgens <- IUCN_allgens[which(IUCN_allgens$number.of.accepted.species > 1 & IUCN_allgens$assessments_per_genus/IUCN_allgens$number.of.accepted.species >= 0.2),]
wilcox.test(with(IUCN_allgens,threats_per_genus_criterionB/assessments_per_genus), with(IUCN_ages_sub, n_threat_criterionB/n_assess), var.equal=F)
mean(with(IUCN_ages_sub, n_threat_criterionB/n_assess))
sd(with(IUCN_ages_sub, n_threat_criterionB/n_assess))/sqrt(length(with(IUCN_ages_sub, n_threat_criterionB/n_assess)))
mean(with(IUCN_allgens,threats_per_genus_criterionB/assessments_per_genus))
sd(with(IUCN_allgens,threats_per_genus_criterionB/assessments_per_genus))/sqrt(length(with(IUCN_allgens,threats_per_genus_criterionB/assessments_per_genus)))

# is there a bias in sampling rare/small-ranged species?
fam_stats <- as.data.frame(cbind( with(IUCN_ages_sub, tapply(n_threat_criterionB, family, sum)/tapply(n_assess, family, sum)),
                                  with(IUCN_ages_sub, tapply(prop.genera.sampled, family, unique)) ))
colnames(fam_stats) <- c('prop.threat','prop.sampled')
PhytoPhylo_familysub <- drop.tip(IUCN_PhytoPhylo_merged$phy,as.character(IUCN_ages_sub[duplicated(IUCN_ages_sub$family),'tip']))
PhytoPhylo_familysub$tip.label <- IUCN_ages_sub[match(PhytoPhylo_familysub$tip.label,IUCN_ages_sub$tip),'family']
fam_stats <- fam_stats[match(PhytoPhylo_familysub$tip.label, rownames(fam_stats)),]
fam_stats$tip <- rownames(fam_stats)
fam_merged <- comparative.data(phy=PhytoPhylo_familysub,data=fam_stats,names.col=tip,vcv=T)
summary(pgls(prop.threat ~ prop.sampled, data=fam_merged,lambda='ML'))





################################################################################
###### 5) HOW SIMILAR ARE GENUS AND SPECIES LEVEL AGES?
################################################################################
# calculate tip ages
PhytoPhylo_nooutliers <- drop.tip(PhytoPhylo_nooutliers,which(duplicated(PhytoPhylo_nooutliers$tip.label)))
ages1_spp <- calculate_tip_ages(PhytoPhylo_nooutliers)
ages1_spp[,2] <- as.numeric(as.character(ages1_spp[,2]))
ages1_spp$genus <- sapply(as.character(ages1_spp$tip),function(x) strsplit(x,split='_')[[1]][[1]])

# find genus-level sampling coverage in genus tree
ages1_spp$number.of.species <- plant_lookup_table$number.of.accepted.species[match(ages1_spp$genus, plant_lookup_table$genus)]
ages1_spp$family <- plant_lookup_table$family[match(ages1_spp$genus, plant_lookup_table$genus)]
ages1_spp$group <- plant_lookup_table$group[match(ages1_spp$genus, plant_lookup_table$genus)]
ages1_spp <- ages1_spp[which(is.na(ages1_spp$number.of.species)==F & ages1_spp$number.of.species > 0),]
ages1_spp$prop.spp.sampled <- as.numeric(table(ages1_spp$genus)[match(ages1_spp$genus, names(table(ages1_spp$genus)))])/ages1_spp$number.of.species
ages1_spp$number.of.genera <- as.numeric(genera_in_families[match(ages1_spp$family, names(genera_in_families))])
ages1_spp$prop.genera.sampled <- as.numeric(genera_per_family_phylo[match(ages1_spp$family, names(genera_per_family_phylo))])/ages1_spp$number.of.genera

# drop species from genera with poor sampling coverage
ages1_spp_sub <- ages1_spp[which(ages1_spp$prop.spp.sampled >= .6 & (ages1_spp$number.of.species > 1 | ages1_spp$number.of.species == 1 & ages1_spp$prop.genera.sampled >= 0.6)),]

# merge with genus ages
ages1_spp_sub_genus <- merge(ages1_spp_sub,ages1,by.x='genus',by.y='tip')
with(ages1_spp_sub_genus, cor.test(log(tip.age.x), log(tip.age.y)))




################################################################################
###### 6) LOAD DATA FOR SPECIES-LEVEL ANALYSIS
################################################################################
#read conifer phylogeny
con_tree <- read.nexus('Leslie et al_AJB 2018.tree')

# read the palm phylogeny
palm_tree <- read.tree('palm_tree.tre')
palm_tree_post <- read.tree('palm_posterior.trees')

# calculate species ages in conifers, merge with IUCN data and clean up names with TPL
species_ages<-calculate_tip_ages(con_tree)
con_genus<-sapply(as.character(species_ages$tip),function(x) strsplit(x,split='_')[[1]][[1]])
con_species<-sapply(as.character(species_ages$tip),function(x) { qq1 <- strsplit(x,split='_')[[1]]; if(length(qq1==2)){return(qq1[2])}; if(length(qq1==1)){return("")}})
con_TPL <- TPL(genus=con_genus,species=con_species)
con_TPL$Merged <- with(con_TPL, paste(New.Genus,New.Species,sep='_'))
con_TPL_IUCN <- merge(con_TPL, IUCN_TPL_riskunique_drop[,c('Merged','Red.List.status','Red.List.criteria','Red.List.criteria.version','Year.assessed')], by='Merged', all.x=T)
con_TPL_IUCN$old_species <- with(con_TPL_IUCN, paste(Genus,Species,sep='_'))
species_ages$old_species <- paste(con_genus,con_species,sep='_')
species_ages_est <- merge(species_ages, con_TPL_IUCN, by = "old_species", all.x=T)
species_ages_est$tip <- as.character(species_ages_est$tip)
species_ages_est$tip.age <- as.numeric(as.character(species_ages_est$tip.age))

# calculate species ages in palms and clean up names with TPL
palm_ages_o<-calculate_tip_ages(palm_tree)
palm_ages_o_all<-lapply(palm_tree_post, calculate_tip_ages)
palm_genera<-unname(sapply(as.character(palm_ages_o$tip),function(x)strsplit(x,split='_')[[1]][1]))
palm_species<-unname(sapply(as.character(palm_ages_o$tip),function(x)strsplit(x,split='_')[[1]][2]))
palm_TPL<-TPL(genus=palm_genera,species=palm_species)
palm_TPL[,'Merged']<-paste(palm_TPL$New.Genus,palm_TPL$New.Species, sep = '_')
palm_TPL_IUCN <- merge(palm_TPL, IUCN_TPL_riskunique_drop[,c('Merged','Red.List.status','Red.List.criteria.version','Red.List.criteria','Year.assessed')], by='Merged', all.x=T)
palm_TPL_IUCN$old_species <- with(palm_TPL_IUCN, paste(Genus,Species,sep='_'))
palm_ages <- merge(palm_ages_o, palm_TPL_IUCN, by.x = "tip", by.y = "old_species", all.x=T)
palm_ages$tip <- as.character(palm_ages$tip)
palm_ages$tip.age <- as.numeric(as.character(palm_ages$tip.age))
palm_ages_all <- lapply(palm_ages_o_all, function(x) {
                        merged_palms <- merge(x, palm_TPL_IUCN, by.x = "tip", by.y = "old_species", all.x=T)
                        merged_palms$tip <- as.character(merged_palms$tip)
                        merged_palms$tip.age <- as.numeric(as.character(merged_palms$tip.age))
                        return(merged_palms)
                        })

#remove synonyms
species_ages_est_a <- species_ages_est[which(species_ages_est$Taxonomic.status=='Accepted'),]
palm_ages <- palm_ages[which(palm_ages$Taxonomic.status=='Accepted'),]
palm_dups <- palm_ages[duplicated(palm_ages$Merged),]$Merged
palm_ages <- palm_ages[-which(palm_ages$Merged %in% palm_dups & (with(palm_ages,paste(Genus,Species,sep="_")) %in% palm_ages$Merged == F)),]
palm_ages_all <- lapply(palm_ages_all, function(x) {
                        tmp1 <- x[which(x$Taxonomic.status=='Accepted'),]
                        palm_dups <- tmp1[duplicated(tmp1$Merged),]$Merged
                        tmp1 <- tmp1[-which(tmp1$Merged %in% palm_dups & (with(tmp1,paste(Genus,Species,sep="_")) %in% tmp1$Merged == F)),]
                        return(tmp1)
                        })

# simplify risk status
species_ages_est_a <- species_ages_est_a[species_ages_est_a$Red.List.status %in% c('EW','EX','DD') == F,]
#species_ages_est_a$Red.List.status <- as.numeric(species_ages_est_a$Red.List.status=="Threatened")
levels(palm_ages$Red.List.status)[levels(palm_ages$Red.List.status)%in%c("CR","EN","VU")] <- "Threatened"
levels(palm_ages$Red.List.status)[levels(palm_ages$Red.List.status)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
palm_ages <- palm_ages[palm_ages$Red.List.status %in% c('EW','EX','DD') == F,]
palm_ages$Red.List.status <- as.numeric(palm_ages$Red.List.status=="Threatened")
palm_ages_all <- lapply(palm_ages_all, function(x) {
                        levels(x$Red.List.status)[levels(x$Red.List.status)%in%c("CR","EN","VU")] <- "Threatened"
                        levels(x$Red.List.status)[levels(x$Red.List.status)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
                        tmp1 <- x[x$Red.List.status %in% c('EW','EX','DD') == F,]
                        tmp1$Red.List.status <- as.numeric(tmp1$Red.List.status=="Threatened")
                        return(tmp1)
                        })





################################################################################
###### 7) TEST ASSOCIATION AMONG AGE, THREAT STATUS, AND RANGE IN CONIFERS
################################################################################
# find tips with >=90% support
keepc <- c('Podocarpus_decumbens','Podocarpus_spinulosus','Podocarpus_glaucus','Podocarpus_latifolius','Podocarpus_milanjianus','Podocarpus_elongatus','Podocarpus_henkelii',
  'Retrophyllum_comptonii','Retrophyllum_minus','Dacrydium_xanthandrum','Dacrydium_magnum','Dacrydium_cupressinum','Dacrycarpus_kinabaluensis','Dacrycarpus_cinctus','Dacrycarpus_cumingii',
  'Acmopyle_pancheri','Acmopyle_sahniana','Pherosphaera_fitzgeraldii','Pherosphaera_hookeriana','Microcachrys_tetragona','Saxegothaea_conspicua','Manoao_colensoi','Lagarostrobos_franklinii','Parasitaxus_usta',
  'Halocarpus_kirkii','Phyllocladus_trichomanoides','Lepidothamnus_laxifolius','Araucaria_muelleri','Araucaria_cunninghamii','Araucaria_hunsteinii','Agathis_australis','Wollemia_nobilis',
  'Juniperus_flaccida','Juniperus_monticola','Juniperus_gracilior','Juniperus_bermudiana','Juniperus_excelsa','Juniperus_polycarpos','Juniperus_thurifera',
  'Hesperocyparis_benthamii','Hesperocyparis_bakeri','Cupressus_funebris','Calocedrus_macrolepis','Calocedrus_decurrens','Platycladus_orientalis','Microbiota_decussata',
  'Chamaecyparis_pisifera','Chamaecyparis_formosensis','Chamaecyparis_obtusa','Chamaecyparis_lawsoniana','Fokienia_hodginsii','Thuja_standishii','Thujopsis_dolabrata','Callitris_endlicheri','Callitris_rhomboidea',
  'Callitris_canescens','Actinostrobus_pyramidalis','Actinostrobus_arenarius','Actinostrobus_acuminatus','Callitris_drummondii','Callitris_columellaris',
  'Callitris_sulcata','Neocallitropsis_pancheri','Widdringtonia_schwarzii','Widdringtonia_cedarbergensis','Widdringtonia_nodiflora','Widdringtonia_whytei','Diselma_archeri','Fitzroya_cupressoides',
  'Libocedrus_bidwillii','Pilgerodendron_uviferum','Austrocedrus_chilensis','Papuacedrus_papuana','Taxodium_distichum','Taxodium_mucronatum','Glyptostrobus_pensilis',
  'Cryptomeria_japonica','Sequoia_sempervirens','Sequoiadendron_giganteum','Metasequoia_glyptostroboides','Athrotaxis_cupressoides','Taiwania_cryptomerioides',
  'Cunninghamia_lanceolata','Cunninghamia_konishii','Taxus_floridana','Pseudotaxus_chienii','Austrotaxus_spicata','Torreya_nucifera','Torreya_taxifolia','Amentotaxus_argotaenia','Amentotaxus_formosana',
  'Cephalotaxus_harringtonii','Cephalotaxus_oliveri','Sciadopitys_verticillata','Pinus_clausa','Pinus_virginiana','Pinus_contorta','Pinus_mugo','Pinus_uncinata','Pinus_merkusii',
  'Pinus_latteri','Pinus_brutia','Pinus_halepensis','Pinus_monticola','Pinus_flexilis','Pinus_ayacahuite','Pinus_strobiformis','Pinus_chiapensis','Pinus_strobus','Pinus_squamata',
  'Pinus_maximartinezii','Pinus_pinceana','Pinus_balfouriana','Picea_engelmannii','Picea_breweriana','Abies_sibirica','Abies_alba',
  'Keteleeria_evelyniana','Nothotsuga_longibracteata','Pseudolarix_amabilis','Cedrus_deodara',
  'Podocarpus_elatus','Podocarpus_polystachyus','Podocarpus_drouynianus','Podocarpus_angustifolius','Podocarpus_ekmanii','Podocarpus_totara','Podocarpus_acutifolius','Agathis_vitiensis',
  'Nageia_formosensis','Nageia_nagi','Retrophyllum_vitiense','Juniperus_saltillensis','Juniperus_durangensis','Juniperus_oxycedrus','Cupressus_austrotibetica','Cupressus_cashmeriana','Cupressus_torulosa',
  'Callitris_preissii','Callitris_verrucosa','Pinus_banksiana','Pinus_peuce','Pinus_monophylla','Pinus_quadrifolia','Pinus_culminicola','Pinus_discolor','Larix_laricina','Larix_kaempferi','Abies_magnifica',
  'Abies_mariesii',
  'Tetraclinis_articulata','Athrotaxis_selaginoides','Athrotaxis_laxifolia',
  # retain outgroups in tree
  'Encephalartos','Zamia','Cycas'
)

# load distribution data
conifer_ranges <- read.csv('coniferSDM566.csv', header = T)
conifer_distru <- read.table('coniferSDM566_GBIFlatitude.txt', header = T)
conifer_distru <- read.csv('latRangeConifers.csv', header = T)
species_data_sub <- merge(species_ages_est_a,conifer_ranges,by.x='Merged',by.y='sp')
species_data_sub <- merge(species_data_sub,conifer_distru,by.x='Merged',by.y='species')
species_data_sub$gbif_recs <- species_data_sub$d.tpos + species_data_sub$d.fneg
species_data_sub$Red.List.status <- species_data_sub$Red.List.status[drop=T]
species_data_sub$accuracy <- with(species_data_sub, (d.tpos+d.tneg)/(d.fpos+d.fneg+d.tpos+d.tneg))

con_data_merged_all <- comparative.data(phy=con_tree,data=species_data_sub[which(species_data_sub$tip %in% keepc == T),c('tip','tip.age','gbif_recs','global.count','Genus','EA_0.25_cells','Red.List.status','Red.List.criteria')],names.col=tip,vcv=T)

# fit individual models
summary( pgls(log(global.count)~scale(log(tip.age)), con_data_merged_all, lambda='ML') )
summary(phyloglm(as.numeric( Red.List.status %in% c('EN','CR','VU') )~scale(log(global.count)) + scale(log(tip.age)),phy=con_data_merged_all$phy,data=con_data_merged_all$data,method="logistic_MPLE"))

# use phylogenetic path analysis
con_data_merged_all$data$status <- as.factor(con_data_merged_all$data$Red.List.status %in% c('EN','CR','VU') & ( grepl('B',con_data_merged_all$data$Red.List.criteria) |
                                                                                                                 grepl("D2",con_data_merged_all$data$Red.List.criteria)|
                                                                                                                 con_data_merged_all$data$Red.List.criteria == 'D1+2')  )
con_data_merged_all$data$s_age <- as.numeric(scale(log(con_data_merged_all$data$tip.age)))
con_data_merged_all$data$s_range <- as.numeric(scale(log(con_data_merged_all$data$global.count)))

con_path_models <- define_model_set(
                                    one = c(status ~ s_age + s_range),
                                    two = c(status ~ s_range, s_range ~ s_age),
                                    three = c(status ~ s_range)
                                    )
con_path_result <- phylo_path(con_path_models, data = con_data_merged_all$data, tree = con_data_merged_all$phy, model = 'lambda', method = 'logistic_MPLE')
R2.lik(con_path_result$d_sep$two$model[[1]])
R2.pred(con_path_result$d_sep$one$model[[1]],phy=con_data_merged_all$phy)





################################################################################
###### 8) SISTER SPECIES ANALYSIS IN CONIFERS
################################################################################
# find sisters in full tree and then use those with contrasting extinction status
sister_extract <- cbind(con_tree$tip,as.character(unlist(sapply(con_tree$tip,getSisters,tree=con_tree,mode='label'))))
sister_extract <- sister_extract[!duplicated(t(apply(sister_extract, 1, sort))), ]
sister_extract <- sister_extract[apply(sister_extract,1,function(x){sum(is.na(as.numeric(x)))})==2,]
sister_extract <- cbind(sister_extract, as.data.frame(t(apply(sister_extract, 1, function(x){  species_ages_est_a[match(x, species_ages_est_a$tip),]$Red.List.status  }))))
sister_extract <- cbind(sister_extract, as.data.frame(t(apply(sister_extract, 1, function(x){  species_ages_est_a[match(x[1:2], species_ages_est_a$tip),]$Red.List.criteria  }))))
colnames(sister_extract) <- c('Sp1','Sp2','Stat1','Stat2','Crit1','Crit2')
sister_extract <- sister_extract[which(is.na(sister_extract$Stat1) == F & is.na(sister_extract$Stat2) == F),]

levels(sister_extract$Stat1)[levels(sister_extract$Stat1)%in%c("CR","EN","VU")] <- "Threatened"
levels(sister_extract$Stat1)[levels(sister_extract$Stat1)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
levels(sister_extract$Stat2)[levels(sister_extract$Stat2)%in%c("CR","EN","VU")] <- "Threatened"
levels(sister_extract$Stat2)[levels(sister_extract$Stat2)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
sister_extract$Stat1 <- as.numeric(sister_extract$Stat1 == 'Threatened' & (grepl('B',sister_extract$Crit1)|grepl("D2",sister_extract$Crit1)|sister_extract$Crit1=='D1+2'))
sister_extract$Stat2 <- as.numeric(sister_extract$Stat2 == 'Threatened' & (grepl('B',sister_extract$Crit2)|grepl("D2",sister_extract$Crit2)|sister_extract$Crit2=='D1+2'))

# sisters with >=90% support
sister_extract <- sister_extract[c(which(sister_extract[,1] %in% keepc)),]

# find the identity of sisters with differing extinction statuses                                                          
sister_extract_nomatch <- sister_extract[which(sister_extract$Stat1 != sister_extract$Stat2),]

# extract characteristics for data frame of sisters
sister_extract_dat <- as.data.frame(t( apply(sister_extract_nomatch, 1, sis_chars, dataf = con_data_merged_all$data, sampling = F) ))
colnames(sister_extract_dat) <- c('Sp1','Sp2','Time','diff_PotRange')
sister_extract_dat[,3:ncol(sister_extract_dat)] <- sapply(sister_extract_dat[,3:ncol(sister_extract_dat)], function(x){ as.numeric(as.character(x))})

# find sisters that are both non-threatened and sample from them 
sister_extract_match <- sister_extract[which(sister_extract$Stat1 == sister_extract$Stat2 & sister_extract$Stat2 == F),]
ran_cors <- lapply(1:1000, function(y){
   sister_extract_dat2 <- as.data.frame(t( apply(sister_extract_match[sample(1:nrow(sister_extract_match), nrow(sister_extract_nomatch), replace=F), ], 1, sis_chars, dataf = con_data_merged_all$data, sampling = T) ))
   colnames(sister_extract_dat2) <- c('Sp1','Sp2','Time','diff_PotRange')
   sister_extract_dat2[,3:ncol(sister_extract_dat2)] <- sapply(sister_extract_dat2[,3:ncol(sister_extract_dat2)], function(x){ as.numeric(as.character(x))})
   output <- sapply(which(colnames(sister_extract_dat2)=='diff_PotRange'):ncol(sister_extract_dat2), function(z){ cor(sister_extract_dat2[,z],log(sister_extract_dat2$Time)) })
   return(output)
  })
ran_cors <- do.call("rbind",ran_cors)
tru_cors <- sapply(which(colnames(sister_extract_dat)=='diff_PotRange'):ncol(sister_extract_dat), function(z){ cor(sister_extract_dat[,z],log(sister_extract_dat$Time),use='pairwise.complete.obs') })
names(tru_cors) <- 'diff_PotRange'


# check if NT in sister_extract_nomatch differs from NT in sister_extract_match
w1 <- apply(sister_extract_nomatch, 1, function(x){  con_data_merged_all$data[match(x[1:2][which(x[3:4] == 0)],rownames(con_data_merged_all$data)),'global.count']   })
w2 <- con_data_merged_all$data[match( as.character(unlist(sister_extract_match[1:2])), rownames(con_data_merged_all$data)),'global.count']
leveneTest(log(c(w1,w2)) ~ c(rep('d',length(w1)),rep('s',length(w2))) )
t.test(log(c(w1,w2)) ~ c(rep('d',length(w1)),rep('s',length(w2))),var.equal=T )


# find sampling coverage in the focal genera and exclude those with <90%
gen_acc <- as.data.frame(cbind(unique(species_ages_est_a$Genus),as.numeric(with(species_ages_est_a, table(Genus)))))
gen_acc$number.of.accepted.species <- plant_lookup_table[match(gen_acc[,1],plant_lookup_table$genus),'number.of.accepted.species']
gen_acc$sampling <- as.numeric(as.character(gen_acc[,2]))/gen_acc$number.of.accepted.species*100
gen_acc_sis <- gen_acc[match(unique(sapply(strsplit(as.character(sister_extract$Sp1),'_'),function(x){x[1]})),gen_acc[,1]),]
sister_extract <- sister_extract[c(which(sapply(strsplit(as.character(sister_extract$Sp1),'_'),function(x){x[1]}) %in% gen_acc_sis[which(gen_acc_sis$sampling >= 90),1])),]

# find the identity of sisters with differing extinction statuses
sister_extract_nomatch <- sister_extract[which(sister_extract$Stat1 != sister_extract$Stat2),]

# extract characteristics for data frame of sisters
sister_extract_dat <- as.data.frame(t( apply(sister_extract_nomatch, 1, sis_chars, dataf = con_data_merged_all$data, sampling = F) ))
colnames(sister_extract_dat) <- c('Sp1','Sp2','Time','diff_PotRange')
sister_extract_dat[,3:ncol(sister_extract_dat)] <- sapply(sister_extract_dat[,3:ncol(sister_extract_dat)], function(x){ as.numeric(as.character(x))})

# find sisters that are both non-threatened and sample from them
sister_extract_match <- sister_extract[which(sister_extract$Stat1 == sister_extract$Stat2 & sister_extract$Stat2 == F),]
ran_cors <- lapply(1:100, function(y){
   sister_extract_dat2 <- as.data.frame(t( apply(sister_extract_match[sample(1:nrow(sister_extract_match), nrow(sister_extract_nomatch), replace=F), ], 1, sis_chars, dataf = con_data_merged_all$data, sampling = T) ))
   colnames(sister_extract_dat2) <- c('Sp1','Sp2','Time','diff_PotRange')
   sister_extract_dat2[,3:ncol(sister_extract_dat2)] <- sapply(sister_extract_dat2[,3:ncol(sister_extract_dat2)], function(x){ as.numeric(as.character(x))})
   output <- sapply(which(colnames(sister_extract_dat2)=='diff_PotRange'):ncol(sister_extract_dat2), function(z){ cor(sister_extract_dat2[,z],log(sister_extract_dat2$Time)) })
   return(output)
  })
ran_cors <- do.call("rbind",ran_cors)
tru_cors <- sapply(which(colnames(sister_extract_dat)=='diff_PotRange'):ncol(sister_extract_dat), function(z){ cor(sister_extract_dat[,z],log(sister_extract_dat$Time),use='pairwise.complete.obs') })
names(tru_cors) <- 'diff_PotRange'





################################################################################
###### 9) RERUN ANALYSIS AT GENUS LEVEL IN CONIFERS
################################################################################
# filter tree
generaPhyndr <- unique(sapply(strsplit(as.character(con_tree$tip.label),'_'),function(x) x[1]))
generaPhyndr.2 <- sapply(generaPhyndr, function(x) x<- paste(x,'_',sep=''))
ii<-sapply(generaPhyndr.2,function(x,y) grep(x,y,fixed=TRUE)[1],y=con_tree$tip.label)
#drop all but one of every genera
Con_GenusTree<-drop.tip(con_tree,setdiff(con_tree$tip.label,con_tree$tip.label[ii]))
#remove full species name and leave only genera as tip label
Con_GenusTree$tip.label<-sapply(strsplit(Con_GenusTree$tip.label,"_"),function(x) x[1])

# assemble dataframe
species_ages[,2] <- as.numeric(as.character(species_ages[,2]))
con_gen_data <- as.data.frame(cbind( names(with(species_ages, tapply(tip.age, sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max))),
                                      as.numeric(with(species_ages, tapply(tip.age, sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max)))
                    ))
con_gen_data[,1] <- as.character(con_gen_data[,1])
con_gen_data[,2] <- as.numeric(as.character(con_gen_data[,2]))
colnames(con_gen_data) <- c("genus","age")
con_gen_data <- merge(con_gen_data,plant_lookup_table,all.x=T,by="genus")

# for each genus calculate diversification index
con_gen_data$r_hat0 <- bd.ms(time=con_gen_data$age, n=con_gen_data$number.of.accepted.species, crown=F, epsilon = 0)
con_gen_data$r_hat50 <- bd.ms(time=con_gen_data$age, n=con_gen_data$number.of.accepted.species, crown=F, epsilon = 0.5)
con_gen_data$r_hat90 <- bd.ms(time=con_gen_data$age, n=con_gen_data$number.of.accepted.species, crown=F, epsilon = 0.9)
con_gen_data_old <- con_gen_data

# merge and simplify risk status
con_TPL_IUCN2 <- con_TPL_IUCN[which(con_TPL_IUCN$Taxonomic.status=="Accepted" & is.na(con_TPL_IUCN$Red.List.status)==F),]
levels(con_TPL_IUCN2$Red.List.status)[levels(con_TPL_IUCN2$Red.List.status)%in%c("CR","EN","VU")] <- "Threatened"
levels(con_TPL_IUCN2$Red.List.status)[levels(con_TPL_IUCN2$Red.List.status)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
con_TPL_IUCN2$Red.List.statusB <- as.numeric(con_TPL_IUCN2$Red.List.status=="Threatened" & (grepl('B',con_TPL_IUCN2$Red.List.criteria)|grepl('D1',con_TPL_IUCN2$Red.List.criteria)|con_TPL_IUCN2$Red.List.criteria == 'D1+2'))

con_TPL_IUCN2$Red.List.status <- as.numeric(con_TPL_IUCN2$Red.List.status=="Threatened")
n_threat <- with(con_TPL_IUCN2, tapply(Red.List.status,Genus[drop=T],sum))
n_threatB <- with(con_TPL_IUCN2, tapply(Red.List.statusB,Genus[drop=T],sum))
n_assess <- with(con_TPL_IUCN2, tapply(Red.List.status,Genus[drop=T],length))

con_gen_data <- con_gen_data[which(con_gen_data$genus %in% names(n_threatB)),]
con_gen_data$n_threat<-as.numeric(n_threatB[match(con_gen_data$genus,names(n_threatB))])
con_gen_data$n_assess<-as.numeric(n_assess[match(con_gen_data$genus,names(n_assess))])

# prepare data for PGLS
con_gen_data_sub <- con_gen_data[which(con_gen_data$r_hat0 != 0 & con_gen_data$n_assess/con_gen_data$number.of.accepted.species >= 0.2 & con_gen_data$number.of.accepted.species > 1),]
Con_GenusTree_IUCNsub <- drop.tip(Con_GenusTree,setdiff(Con_GenusTree$tip.label,con_gen_data_sub$genus))
con_gen_data_sub <- con_gen_data_sub[match(Con_GenusTree_IUCNsub$tip.label, con_gen_data_sub$genus),]

# re-run PGLS weighted by the proportion of species that are assessed
gls_wts <- sqrt(con_gen_data_sub$n_assess/con_gen_data_sub$number.of.accepted.species)
mod1b <- gls(I(n_threat/n_assess) ~ log(r_hat0),  data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod2b <- gls(I(n_threat/n_assess) ~ log(r_hat50), data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod3b <- gls(I(n_threat/n_assess) ~ log(r_hat90), data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))

# prepare data for PGLS
con_gen_data_sub2 <- con_gen_data[which(con_gen_data$n_assess/con_gen_data$number.of.accepted.species >= 0.2 & con_gen_data$number.of.accepted.species > 0),]
Con_GenusTree_IUCNsub2 <- drop.tip(Con_GenusTree,setdiff(Con_GenusTree$tip.label,con_gen_data_sub2$genus))
con_gen_data_sub2 <- con_gen_data_sub2[match(Con_GenusTree_IUCNsub2$tip.label, con_gen_data_sub2$genus),]
gls_wts2 <- sqrt(con_gen_data_sub2$n_assess/con_gen_data_sub2$number.of.accepted.species)
mod4b <- gls(I(n_threat/n_assess) ~ log(age), data=con_gen_data_sub2, correlation=corPagel(0,Con_GenusTree_IUCNsub2,fixed=T),method="ML",weights=varFixed(~1/gls_wts2))





################################################################################
###### 10) CHECK BIAS IN THE CONIFER DATASET
################################################################################
# table S1
summary( pgls(log(gbif_recs)~scale(log(tip.age))*as.numeric( status ), con_data_merged_all, lambda='ML') )

# discussion point
PhytoPhylo_nooutliers_nodup <- drop.tip(PhytoPhylo_nooutliers,which(duplicated(PhytoPhylo_nooutliers$tip.label)))
ages2 <- calculate_tip_ages(PhytoPhylo_nooutliers_nodup)
ages2[,2] <- as.numeric(as.character(ages2[,2]))
leveneTest(log(c(con_data_merged_all$data$tip.age,as.numeric(as.character(ages2$tip.age)))), c(rep("data",dim(con_data_merged_all$data)[1]),rep("all",dim(ages2)[1])))
t.test(log(con_data_merged_all$data$tip.age),log(ages2[,2] ), var.equal=T)

# table S4
# species age distribution across tree
# clean-up original "full" data
con_TPL2 <- con_TPL
con_TPL2$old_species <- with(con_TPL2, paste(Genus,Species,sep='_'))
species_ages_allcons <- merge(species_ages, con_TPL2, by.x = "tip", by.y = "old_species", all.x=T)
tmp_con <- species_ages_allcons[which(species_ages_allcons$Taxonomic.status=='Accepted'),]
leveneTest(log(c(con_data_merged_all$data$tip.age,as.numeric(as.character(tmp_con$tip.age)))), c(rep("data",dim(con_data_merged_all$data)[1]),rep("all",dim(tmp_con)[1])))
t.test(log(con_data_merged_all$data$tip.age),log(as.numeric(as.character(tmp_con$tip.age))), var.equal=F)
mean(log(con_data_merged_all$data$tip.age))
sd(log(con_data_merged_all$data$tip.age))/sqrt(length(con_data_merged_all$data$tip.age))
mean(log(as.numeric(as.character(tmp_con$tip.age))))
sd(log(as.numeric(as.character(tmp_con$tip.age))))/sqrt(length(tmp_con$tip.age))

# genus age distribution
con_gen_data_old <- con_gen_data_old[is.na(con_gen_data_old$number.of.accepted.species)==F,]
con_gen_data_old <- con_gen_data_old[which(con_gen_data_old$r_hat0 != 0 & con_gen_data_old$number.of.accepted.species > 1),]
leveneTest(log(c(con_gen_data_old$age,as.numeric(as.character(con_gen_data_sub$age)))), c(rep("data",dim(con_gen_data_old)[1]),rep("all",dim(con_gen_data_sub)[1])))
t.test(log(con_gen_data_sub$age),log(as.numeric(as.character(con_gen_data_old$age))), var.equal=T)
mean(log(con_gen_data_sub$age))
sd(log(con_gen_data_sub$age))/sqrt(length(con_gen_data_sub$age))
mean(log(as.numeric(as.character(con_gen_data_old$age))))
sd(log(as.numeric(as.character(con_gen_data_old$age))))/sqrt(length(con_gen_data_old$age))

# species per genus
leveneTest(log(c(con_gen_data_old$number.of.accepted.species,as.numeric(as.character(con_gen_data_sub$number.of.accepted.species)))), c(rep("data",dim(con_gen_data_old)[1]),rep("all",dim(con_gen_data_sub)[1])))
t.test(log(con_gen_data_sub$number.of.accepted.species),log(as.numeric(as.character(con_gen_data_old$number.of.accepted.species))), var.equal=T)
mean(log(con_gen_data_sub$number.of.accepted.species))
sd(log(con_gen_data_sub$number.of.accepted.species))/sqrt(length(con_gen_data_sub$number.of.accepted.species))
mean(log(as.numeric(as.character(con_gen_data_old$number.of.accepted.species))))
sd(log(as.numeric(as.character(con_gen_data_old$number.of.accepted.species))))/sqrt(length(con_gen_data_old$number.of.accepted.species))

# species per genus
leveneTest(log(c(con_gen_data_old$r_hat50,as.numeric(as.character(con_gen_data_sub$r_hat50)))), c(rep("data",dim(con_gen_data_old)[1]),rep("all",dim(con_gen_data_sub)[1])))
t.test(log(con_gen_data_sub$r_hat50),log(as.numeric(as.character(con_gen_data_old$r_hat50))), var.equal=T)
mean(log(con_gen_data_sub$r_hat50))
sd(log(con_gen_data_sub$r_hat50))/sqrt(length(con_gen_data_sub$r_hat50))
mean(log(as.numeric(as.character(con_gen_data_old$r_hat50))))
sd(log(as.numeric(as.character(con_gen_data_old$r_hat50))))/sqrt(length(con_gen_data_old$r_hat50))





################################################################################
###### 11) TEST ASSOCIATION AMONG AGE, THREAT STATUS, AND RANGE IN PALMS
################################################################################
# load TTR and GBIF raw data
palm_stats <- read.csv('palm_stats.csv')
palm_ages_range <- merge(palm_ages,palm_stats,by.x='Merged',by.y='sp')
palm_ages_range$gbif_recs <- with(palm_ages_range,d.tpos + d.fneg)

# drop those without threat data and merge with ranges
palms_IUCN_drop <- palm_ages[is.na(palm_ages$Red.List.status)==F,]
palms_IUCN_drop_range <- merge(palms_IUCN_drop,palm_ages_range[,c('Merged','EA_0.25_cells','global.count','gbif_recs')],by='Merged')
palms_dropped_tree <- drop.tip(palm_tree, setdiff(palm_tree$tip, palms_IUCN_drop_range$tip) )
palms_dropped_tree$tip.label <- palms_IUCN_drop_range[match(palms_dropped_tree$tip.label, palms_IUCN_drop_range$tip),'Merged']
palms_IUCN_drop_range <- palms_IUCN_drop_range[match(palms_dropped_tree$tip.label,palms_IUCN_drop_range$Merged),]
rownames(palms_IUCN_drop_range) <- palms_IUCN_drop_range$Merged

palm_data_merged <- comparative.data(phy=palm_tree,data=palms_IUCN_drop_range[,c('tip','tip.age','Red.List.status','Red.List.criteria','gbif_recs','global.count')],names.col=tip,vcv=T)

# use phylogenetic path analysis
palm_path_models <- define_model_set(
                                    one = c(Red.List.status ~ s_age + s_range),
                                    two = c(Red.List.status ~ s_range, s_range ~ s_age),
                                    three = c(Red.List.status ~ s_range)
                                    )

# run analysis across posterior of trees
palm_path_result_all <- sapply(1:length(palm_tree_post), function(x){
               tmp1_data <- merge(palm_ages_all[[x]],palm_stats,by.x='Merged',by.y='sp')
               tmp1_data$gbif_recs <- with(tmp1_data,d.tpos + d.fneg)
               tmp1_data <- tmp1_data[is.na(tmp1_data$Red.List.status)==F,]
               tmp1_data <- merge(palms_IUCN_drop,tmp1_data[,c('Merged','EA_0.25_cells','global.count','gbif_recs')],by='Merged')
               tmp1_tree <- drop.tip(palm_tree_post[[x]], setdiff(palm_tree_post[[x]]$tip, tmp1_data$tip) )
               tmp1_tree$tip.label <- tmp1_data[match(tmp1_tree$tip.label, tmp1_data$tip),'Merged']
               tmp1_data <- tmp1_data[match(tmp1_tree$tip.label,tmp1_data$Merged),]
               rownames(tmp1_data) <- tmp1_data$Merged
               tmp1_data$s_age <- as.numeric(scale(log(palm_data_merged$data$tip.age)))
               tmp1_data$s_range <- as.numeric(scale(log(palm_data_merged$data$global.count)))
               tmp1_data$Red.List.status <- as.factor(palm_data_merged$data$Red.List.status==1 & ( grepl("B", palm_data_merged$data$Red.List.criteria) |
                                                                                                   grepl("D2", palm_data_merged$data$Red.List.criteria) | palm_data_merged$data$Red.List.criteria == 'D1+2'))
               ppr <- phylo_path2(palm_path_models, data = tmp1_data, tree = tmp1_tree, model = 'lambda', method = 'logistic_MPLE')

               return( average(ppr) )

               })




# check sampling effect for table S1
gbif_rec_mod <- sapply(1:length(palm_tree_post), function(x){
               tmp1_data <- merge(palm_ages_all[[x]],palm_stats,by.x='Merged',by.y='sp')
               tmp1_data$gbif_recs <- with(tmp1_data,d.tpos + d.fneg)
               tmp1_data <- tmp1_data[is.na(tmp1_data$Red.List.status)==F,]
               tmp1_data <- merge(palms_IUCN_drop,tmp1_data[,c('Merged','EA_0.25_cells','global.count','gbif_recs')],by='Merged')
               tmp1_tree <- drop.tip(palm_tree_post[[x]], setdiff(palm_tree_post[[x]]$tip, tmp1_data$tip) )
               tmp1_tree$tip.label <- tmp1_data[match(tmp1_tree$tip.label, tmp1_data$tip),'Merged']
               tmp1_data <- tmp1_data[match(tmp1_tree$tip.label,tmp1_data$Merged),]
               rownames(tmp1_data) <- tmp1_data$Merged
               tmp1_data$s_age <- as.numeric(scale(log(palm_data_merged$data$tip.age)))
               tmp1_data$s_range <- as.numeric(scale(log(palm_data_merged$data$global.count)))
               tmp1_data$Red.List.status <- as.factor(palm_data_merged$data$Red.List.status==1 & (grepl("B", palm_data_merged$data$Red.List.criteria) |
                                                                                                  grepl("D2", palm_data_merged$data$Red.List.criteria) | palm_data_merged$data$Red.List.criteria == 'D1+2') )
               tmp1_merged <- comparative.data(phy=tmp1_tree,data=tmp1_data[,c('Merged','gbif_recs','Red.List.status','tip.age')],names.col=Merged,vcv=T)
               pgls1 <- pgls(log(gbif_recs)~scale(log(tip.age))*Red.List.status , tmp1_merged, lambda='ML')
               return( as.vector(summary(pgls1)$coefficients[-1,1:2]) )
               })
mu_hat_age  <- sum((1/(455*1000)) * 455*sapply(seq(1,length(gbif_rec_mod),by=6),function(x){gbif_rec_mod[[x]]}))
sd_pool_age <-  sqrt( sum(sapply(seq(4,length(gbif_rec_mod),by=6),function(x){gbif_rec_mod[[x]]})^2)/(length(gbif_rec_mod)/6) )
mu_hat_threat  <- sum((1/(455*1000)) * 455*sapply(seq(2,length(gbif_rec_mod),by=6),function(x){gbif_rec_mod[[x]]}))
sd_pool_threat <-  sqrt( sum(sapply(seq(5,length(gbif_rec_mod),by=6),function(x){gbif_rec_mod[[x]]})^2)/(length(gbif_rec_mod)/6) )
mu_hat_threat_age  <- sum((1/(455*1000)) * 455*sapply(seq(3,length(gbif_rec_mod),by=6),function(x){gbif_rec_mod[[x]]}))
sd_pool_threat_age <-  sqrt( sum(sapply(seq(6,length(gbif_rec_mod),by=6),function(x){gbif_rec_mod[[x]]})^2)/(length(gbif_rec_mod)/6) )
2*(1 - pt(abs(mu_hat_age/sd_pool_age), 455))
2*(1 - pt(abs(mu_hat_threat/sd_pool_threat), 455))
2*(1 - pt(abs(mu_hat_threat_age/sd_pool_threat_age), 455))

# get range of observed values
status_e_fnx <- function(x){
                            srange_est <- con_path_result$d_sep$three$model[[1]]$coef['(Intercept)'] +
                                          average(con_path_result)$coef['s_age','s_range']*(log(x)-mean(log(con_data_merged_all$data$tip.age)))/sd(log(con_data_merged_all$data$tip.age))
                            status_est <- plogis( con_path_result$d_sep$three$model[[2]]$coef['(Intercept)'] +
                                                  average(con_path_result)$coef['s_range','status']*srange_est )
                            return(as.numeric(status_est))
                           }
sapply(range(con_data_merged_all$data$tip.age),status_e_fnx)





################################################################################
###### 12) SISTER SPECIES ANALYSIS IN PALMS
################################################################################
# run across entire posterior
ssa_p <- sapply(1:length(palm_tree_post), function(z){

      # find sisters with contrasting extinction status
      sister_extract_p <- cbind(palm_data_merged$phy$tip,as.character(unlist(sapply(palm_data_merged$phy$tip,getSisters,tree=palm_tree_post[[z]],mode='label'))))
      sister_extract_p <- sister_extract_p[!duplicated(t(apply(sister_extract_p, 1, sort))), ]
      sister_extract_p <- sister_extract_p[apply(sister_extract_p,1,function(x){sum(is.na(as.numeric(x)))})==2,]

      # merge range data with risk status
      tmp1_data <- merge(palm_ages_all[[z]],palm_stats,by.x='Merged',by.y='sp')
      tmp1_data <- tmp1_data[is.na(tmp1_data$Red.List.status)==F,]

      sister_extract_p <- cbind(sister_extract_p, as.data.frame(t(apply(sister_extract_p, 1, function(x){  tmp1_data[match(x, with(tmp1_data,paste(Genus,Species,sep="_"))),]$Red.List.status  }))))
      sister_extract_p <- cbind(sister_extract_p, as.data.frame(t(apply(sister_extract_p, 1, function(x){  tmp1_data[match(x[1:2], with(tmp1_data,paste(Genus,Species,sep="_"))),]$Red.List.criteria  }))))
      colnames(sister_extract_p) <- c('Sp1','Sp2','Stat1','Stat2','Crit1','Crit2')
      sister_extract_p$Stat1 <- as.numeric(sister_extract_p$Stat1 == 1 & (grepl('B',sister_extract_p$Crit1)|grepl("D2",sister_extract_p$Crit1)|sister_extract_p$Crit1=='D1+2'))
      sister_extract_p$Stat2 <- as.numeric(sister_extract_p$Stat2 == 1 & (grepl('B',sister_extract_p$Crit2)|grepl("D2",sister_extract_p$Crit2)|sister_extract_p$Crit2=='D1+2'))

      sister_extract_p <- sister_extract_p[which(is.na(sister_extract_p$Stat1) == F & is.na(sister_extract_p$Stat2) == F),]
      sister_extract_p_nomatch <- sister_extract_p[which(sister_extract_p$Stat1 != sister_extract_p$Stat2),]
      combos <- with(sister_extract_p, table(Stat1,Stat2))

      # extract characteristics for data frame of sisters
      rownames(tmp1_data) <- tmp1_data$tip
      sister_extract_p_nomatch[,1] <- as.character(sister_extract_p_nomatch[,1])
      sister_extract_p_nomatch[,2] <- as.character(sister_extract_p_nomatch[,2])
      sister_extract_dat_p <- as.data.frame(t( sapply(1:nrow(sister_extract_p_nomatch), function(x){sis_chars(sister_extract_p_nomatch[x,], dataf = tmp1_data, sampling = F)}) ))
      colnames(sister_extract_dat_p) <- c('Sp1','Sp2','Time','diff_PotRange')
      sister_extract_dat_p[,3:ncol(sister_extract_dat_p)] <- sapply(sister_extract_dat_p[,3:ncol(sister_extract_dat_p)], function(x){ as.numeric(as.character(x))})

      # find sisters that are both non-threatened and sample from them
      sister_extract_match_p <- sister_extract_p[which(sister_extract_p$Stat1 == sister_extract_p$Stat2 & sister_extract_p$Stat2 == 0),]
      sister_extract_match_p[,1] <- as.character(sister_extract_match_p[,1])
      sister_extract_match_p[,2] <- as.character(sister_extract_match_p[,2])
      sister_extract_dat2 <- as.data.frame(t( sapply(1:nrow(sister_extract_match_p), function(x){sis_chars(sister_extract_match_p[x,], dataf = tmp1_data, sampling = T)}) ))
      colnames(sister_extract_dat2) <- c('Sp1','Sp2','Time','diff_PotRange')
      sister_extract_dat2[,3:ncol(sister_extract_dat2)] <- sapply(sister_extract_dat2[,3:ncol(sister_extract_dat2)], function(x){ as.numeric(as.character(x))})
      ran_cors_p <- sapply(which(colnames(sister_extract_dat2)=='diff_PotRange'):ncol(sister_extract_dat2), function(z){ cor(sister_extract_dat2[,z],log(sister_extract_dat2$Time)) })
      tru_cors_p <- sapply(which(colnames(sister_extract_dat_p)=='diff_PotRange'):ncol(sister_extract_dat_p), function(z){ cor(sister_extract_dat_p[,z],log(sister_extract_dat_p$Time)) })
      names(tru_cors_p) <- c('diff_PotRange')

      return( list(ran_cors_p, tru_cors_p, combos) )
      })
1-sum(sapply(seq(1,3000,3), function(x) { abs(ssa_p[[x]])>abs(ssa_p[[x+1]]) } ))/1000
ssa_p_combos <- sapply(seq(3,3000,3), function(x) { (ssa_p[[x]][1,2]+ssa_p[[x]][2,1]) / ssa_p[[x]][2,2] } )
quantile(ssa_p_combos, prob=c(0.5,0.025,0.975))


# check if T in sister_extract_nomatch differs from T in sister_extract_match
ssa_rang <- sapply(1:length(palm_tree_post), function(z){

      # find sisters with contrasting extinction status
      sister_extract_p <- cbind(palm_data_merged$phy$tip,as.character(unlist(sapply(palm_data_merged$phy$tip,getSisters,tree=palm_tree_post[[z]],mode='label'))))
      sister_extract_p <- sister_extract_p[!duplicated(t(apply(sister_extract_p, 1, sort))), ]
      sister_extract_p <- sister_extract_p[apply(sister_extract_p,1,function(x){sum(is.na(as.numeric(x)))})==2,]

      # merge range data with risk status
      tmp1_data <- merge(palm_ages_all[[z]],palm_stats,by.x='Merged',by.y='sp')
      tmp1_data <- tmp1_data[is.na(tmp1_data$Red.List.status)==F,]

      sister_extract_p <- cbind(sister_extract_p, as.data.frame(t(apply(sister_extract_p, 1, function(x){  tmp1_data[match(x, with(tmp1_data,paste(Genus,Species,sep="_"))),]$Red.List.status  }))))
      sister_extract_p <- cbind(sister_extract_p, as.data.frame(t(apply(sister_extract_p, 1, function(x){  tmp1_data[match(x[1:2], with(tmp1_data,paste(Genus,Species,sep="_"))),]$Red.List.criteria  }))))
      colnames(sister_extract_p) <- c('Sp1','Sp2','Stat1','Stat2','Crit1','Crit2')
      sister_extract_p$Stat1 <- as.numeric(sister_extract_p$Stat1 == 1 & (grepl('B',sister_extract_p$Crit1)|grepl("D2",sister_extract_p$Crit1)|sister_extract_p$Crit1=='D1+2'))
      sister_extract_p$Stat2 <- as.numeric(sister_extract_p$Stat2 == 1 & (grepl('B',sister_extract_p$Crit2)|grepl("D2",sister_extract_p$Crit2)|sister_extract_p$Crit2=='D1+2'))
      sister_extract_p <- sister_extract_p[which(is.na(sister_extract_p$Stat1) == F & is.na(sister_extract_p$Stat2) == F),]
      sister_extract_p_nomatch <- sister_extract_p[which(sister_extract_p$Stat1 != sister_extract_p$Stat2),]
      combos <- with(sister_extract_p, table(Stat1,Stat2))

      # extract characteristics for data frame of sisters
      rownames(tmp1_data) <- tmp1_data$tip
      sister_extract_p_nomatch[,1] <- as.character(sister_extract_p_nomatch[,1])
      sister_extract_p_nomatch[,2] <- as.character(sister_extract_p_nomatch[,2])
      sister_extract_dat_p <- as.data.frame(t( sapply(1:nrow(sister_extract_p_nomatch), function(x){sis_chars(sister_extract_p_nomatch[x,], dataf = tmp1_data, sampling = F)}) ))
      colnames(sister_extract_dat_p) <- c('Sp1','Sp2','Time','diff_PotRange')
      sister_extract_dat_p[,3:ncol(sister_extract_dat_p)] <- sapply(sister_extract_dat_p[,3:ncol(sister_extract_dat_p)], function(x){ as.numeric(as.character(x))})

      # find sisters that are both threatened and sample from them - we can't do non-threatened because there are fewer options than contrasting status
      sister_extract_match_p <- sister_extract_p[which(sister_extract_p$Stat1 == sister_extract_p$Stat2 & sister_extract_p$Stat2 == 1),]
      sister_extract_match_p[,1] <- as.character(sister_extract_match_p[,1])
      sister_extract_match_p[,2] <- as.character(sister_extract_match_p[,2])

      w1 <- apply(sister_extract_p_nomatch, 1, function(x){  tmp1_data[match(x[1:2][which(x[3:4] == 0)],rownames(tmp1_data)),'global.count']   })
      w2 <- tmp1_data[match( as.character(unlist(sister_extract_match_p[1:2])), rownames(tmp1_data)),'global.count']
      LT_p <- leveneTest(log(c(w1,w2)) ~ as.factor(c(rep('d',length(w1)),rep('s',length(w2)))) )[[3]][[1]]

      if(LT_p > 0.05) {
                        tt_p <- t.test(log(c(w1,w2)) ~ c(rep('d',length(w1)),rep('s',length(w2))),var.equal=T)
                      } else {
                        tt_p <- t.test(log(c(w1,w2)) ~ c(rep('d',length(w1)),rep('s',length(w2))),var.equal=F)
                      }
      return( c(tt_p$statistic, tt_p$parameter, tt_p$p.value) )
      })
apply(ssa_rang, 1, quantile, prob=c(0.5,0.025,0.975))





################################################################################
###### 13) RERUN ANALYSIS AT GENUS LEVEL IN PALMS
################################################################################
palms_IUCN_drop$Red.List.statusB <- as.numeric(palms_IUCN_drop$Red.List.status==1 & (grepl('B',palms_IUCN_drop$Red.List.criteria)|grepl('D1',palms_IUCN_drop$Red.List.criteria)|palms_IUCN_drop$Red.List.criteria == 'D1+2'))

# filter tree
gla_p <- sapply(1:length(palm_tree_post), function(z){
                  palms_dropped_tree <- drop.tip(palm_tree_post[[z]], setdiff(palm_tree_post[[z]]$tip, palms_IUCN_drop_range$tip) )
                  generaPhyndr <- unique(sapply(strsplit(as.character(palms_dropped_tree$tip.label),'_'),function(x) x[1]))
                  generaPhyndr.2 <- sapply(generaPhyndr, function(x) x<- paste(x,'_',sep=''))
                  ii<-sapply(generaPhyndr.2,function(x,y) grep(x,y,fixed=TRUE)[1],y=palms_dropped_tree$tip.label)
                  #drop all but one of every genera
                  Palm_GenusTree<-drop.tip(palms_dropped_tree,setdiff(palms_dropped_tree$tip.label,palms_dropped_tree$tip.label[ii]))
                  #remove full species name and leave only genera as tip label
                  Palm_GenusTree$tip.label<-sapply(strsplit(Palm_GenusTree$tip.label,"_"),function(x) x[1])
                  #assemble dataframe
                  palm_gen_data <- as.data.frame(cbind( names(with(palm_ages_o_all[[z]], tapply(as.numeric(as.character(tip.age)), sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max))),
                                      as.numeric(with(palm_ages_o_all[[z]], tapply(as.numeric(as.character(tip.age)), sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max)))
                    ))
                  palm_gen_data[,1] <- as.character(palm_gen_data[,1])
                  palm_gen_data[,2] <- as.numeric(as.character(palm_gen_data[,2]))
                  colnames(palm_gen_data) <- c("genus","age")
                  palm_gen_data <- merge(palm_gen_data,plant_lookup_table,all.x=T,by="genus")
                  n_threat <- with(palms_IUCN_drop, tapply(Red.List.statusB,Genus[drop=T],sum))
                  n_assess <- with(palms_IUCN_drop, tapply(Red.List.status,Genus[drop=T],length))

                  # derive spp counts from tree rather than TPL (some changes don't appear in TPL)
                  tree_based_spp_counts <- table(sapply(strsplit(as.character(palm_tree_post[[z]]$tip.label),'_'),function(x) x[1]))
                  palm_gen_data$tree_based_spp_counts <- tree_based_spp_counts[match(names(tree_based_spp_counts),palm_gen_data$genus)]

                  # for each genus calculate diversification index
                  palm_gen_data$r_hat0 <- bd.ms(time=palm_gen_data$age, n=palm_gen_data$tree_based_spp_counts, crown=F, epsilon = 0)
                  palm_gen_data$r_hat50 <- bd.ms(time=palm_gen_data$age, n=palm_gen_data$tree_based_spp_counts, crown=F, epsilon = 0.5)
                  palm_gen_data$r_hat90 <- bd.ms(time=palm_gen_data$age, n=palm_gen_data$tree_based_spp_counts, crown=F, epsilon = 0.9)

                  # intersect with IUCN data
                  palm_gen_data_old <- palm_gen_data
                  palm_gen_data <- palm_gen_data[which(palm_gen_data$genus %in% names(n_threat)),]
                  palm_gen_data$n_threat<-as.numeric(n_threat[match(palm_gen_data$genus,names(n_threat))])
                  palm_gen_data$n_assess<-as.numeric(n_assess[match(palm_gen_data$genus,names(n_assess))])

                  #run PGLS
                  palm_gen_data_sub <- palm_gen_data[which(palm_gen_data$r_hat0 != 0 & palm_gen_data$n_assess/palm_gen_data$tree_based_spp_counts >= 0.2 & palm_gen_data$tree_based_spp_counts > 1),]
                  Palm_GenusTree_IUCNsub <- drop.tip(Palm_GenusTree,setdiff(Palm_GenusTree$tip.label,palm_gen_data_sub$genus))
                  palm_gen_data_sub <- palm_gen_data_sub[match(Palm_GenusTree_IUCNsub$tip.label, palm_gen_data_sub$genus),]
                  gls_wts <<- as.numeric(sqrt(palm_gen_data_sub$n_assess/palm_gen_data_sub$tree_based_spp_counts))


                  mod1p <- tryCatch({ gls(I(n_threat/n_assess) ~ log(r_hat0), data=palm_gen_data_sub, correlation=corPagel(0,Palm_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts)) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                  if(is.null(mod1p)){
                    mod1p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(r_hat0), data=palm_gen_data_sub, correlation=corPagel(y,Palm_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  if(all(sapply(mod1p, is.null) == T)){
                      mod1p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(r_hat0), data=palm_gen_data_sub, correlation=corPagel(y,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                    }
                  mod1p <- mod1p[which(sapply(mod1p, is.null) == F)][[which.min(sapply(mod1p[which(sapply(mod1p, is.null) == F)],AIC))]]
                  }
                  mod2p <- tryCatch({ gls(I(n_threat/n_assess) ~ log(r_hat50), data=palm_gen_data_sub, correlation=corPagel(0,Palm_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts)) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                  if(is.null(mod2p)){
                    mod2p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(r_hat50), data=palm_gen_data_sub, correlation=corPagel(y,Palm_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  if(all(sapply(mod2p, is.null) == T)){
                      mod2p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(r_hat50), data=palm_gen_data_sub, correlation=corPagel(y,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  }
                  mod2p <- mod2p[which(sapply(mod2p, is.null) == F)][[which.min(sapply(mod2p[which(sapply(mod2p, is.null) == F)],AIC))]]
                  }
                  mod3p <- tryCatch({ gls(I(n_threat/n_assess) ~ log(r_hat90), data=palm_gen_data_sub, correlation=corPagel(0,Palm_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts)) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                  if(is.null(mod3p)){
                    mod3p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(r_hat90), data=palm_gen_data_sub, correlation=corPagel(y,Palm_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  if(all(sapply(mod3p, is.null) == T)){
                      mod3p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(r_hat90), data=palm_gen_data_sub, correlation=corPagel(y,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  }
                  mod3p <- mod3p[which(sapply(mod3p, is.null) == F)][[which.min(sapply(mod3p[which(sapply(mod3p, is.null) == F)],AIC))]]
                  }
                  palm_gen_data_sub2 <- palm_gen_data[which(palm_gen_data$n_assess/palm_gen_data$tree_based_spp_counts >= 0.2 & palm_gen_data$tree_based_spp_counts > 0),]
                  Palm_GenusTree_IUCNsub2 <- drop.tip(Palm_GenusTree,setdiff(Palm_GenusTree$tip.label,palm_gen_data_sub2$genus))
                  palm_gen_data_sub2 <- palm_gen_data_sub2[match(Palm_GenusTree_IUCNsub2$tip.label, palm_gen_data_sub2$genus),]
                  gls_wts2 <<- as.numeric(sqrt(palm_gen_data_sub2$n_assess/palm_gen_data_sub2$tree_based_spp_counts))
                  mod4p <- tryCatch({ gls(I(n_threat/n_assess) ~ log(age), data=palm_gen_data_sub2, correlation=corPagel(0,Palm_GenusTree_IUCNsub2,fixed=F),method="ML",weights=varFixed(~1/gls_wts2)) }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                  if(is.null(mod4p)){
                    mod4p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(age), data=palm_gen_data_sub2, correlation=corPagel(y,Palm_GenusTree_IUCNsub2,fixed=F),method="ML",weights=varFixed(~1/gls_wts2))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  if(all(sapply(mod4p, is.null) == T)){
                      mod4p <- lapply( seq(0,1,length.out=11), function(y){tryCatch({gls(I(n_threat/n_assess) ~ log(age), data=palm_gen_data_sub2, correlation=corPagel(y,Palm_GenusTree_IUCNsub2,fixed=T),method="ML",weights=varFixed(~1/gls_wts2))},error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) })
                  }
                    mod4p <- mod4p[which(sapply(mod4p, is.null) == F)][[which.min(sapply(mod4p[which(sapply(mod4p, is.null) == F)],AIC))]]
                  }
                  print(z)
                  return(c( summary(mod1p)$tTable[2,c('t-value','p-value')],
                            summary(mod2p)$tTable[2,c('t-value','p-value')],
                            summary(mod3p)$tTable[2,c('t-value','p-value')],
                            summary(mod4p)$tTable[2,c('t-value','p-value')],
                            nrow(palm_gen_data_sub), nrow(palm_gen_data_sub2)) )
                 })
apply(gla_p, 1, median)





################################################################################
###### 14) CHECK BIAS IN THE PALM DATASET
################################################################################
palm_biases <- sapply(1:length(palm_ages_all), function(z){

#prepare data at spp-level
tmp1_data <- merge(palm_ages_all[[z]],palm_stats,by.x='Merged',by.y='sp')
tmp1_data <- tmp1_data[is.na(tmp1_data$Red.List.status)==F,]
tmp1_tree <- drop.tip(palm_tree_post[[z]], setdiff(palm_tree_post[[z]]$tip, tmp1_data$tip) )
tmp1_tree$tip.label <- tmp1_data[match(tmp1_tree$tip.label, tmp1_data$tip),'Merged']
tmp1_data <- tmp1_data[match(tmp1_tree$tip.label,tmp1_data$Merged),]
rownames(tmp1_data) <- tmp1_data$Merged

#prepare data at genus-level
palms_dropped_tree <- drop.tip(palm_tree_post[[z]], setdiff(palm_tree_post[[z]]$tip, palms_IUCN_drop_range$tip) )
generaPhyndr <- unique(sapply(strsplit(as.character(palms_dropped_tree$tip.label),'_'),function(x) x[1]))
generaPhyndr.2 <- sapply(generaPhyndr, function(x) x<- paste(x,'_',sep=''))
ii<-sapply(generaPhyndr.2,function(x,y) grep(x,y,fixed=TRUE)[1],y=palms_dropped_tree$tip.label)
Palm_GenusTree<-drop.tip(palms_dropped_tree,setdiff(palms_dropped_tree$tip.label,palms_dropped_tree$tip.label[ii]))
Palm_GenusTree$tip.label<-sapply(strsplit(Palm_GenusTree$tip.label,"_"),function(x) x[1])
palm_gen_data <- as.data.frame(cbind( names(with(palm_ages_o_all[[z]], tapply(as.numeric(as.character(tip.age)), sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max))),
                                      as.numeric(with(palm_ages_o_all[[z]], tapply(as.numeric(as.character(tip.age)), sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max)))
                    ))
palm_gen_data[,1] <- as.character(palm_gen_data[,1])
palm_gen_data[,2] <- as.numeric(as.character(palm_gen_data[,2]))
colnames(palm_gen_data) <- c("genus","age")
palm_gen_data <- merge(palm_gen_data,plant_lookup_table,all.x=T,by="genus")
n_threat <- with(palms_IUCN_drop, tapply(Red.List.statusB,Genus[drop=T],sum))
n_assess <- with(palms_IUCN_drop, tapply(Red.List.status,Genus[drop=T],length))

# derive spp counts from tree rather than TPL (some changes don't appear in TPL)
tree_based_spp_counts <- table(sapply(strsplit(as.character(palm_tree_post[[z]]$tip.label),'_'),function(x) x[1]))
palm_gen_data$tree_based_spp_counts <- tree_based_spp_counts[match(names(tree_based_spp_counts),palm_gen_data$genus)]

# for each genus calculate diversification index
palm_gen_data$r_hat0 <- bd.ms(time=palm_gen_data$age, n=palm_gen_data$tree_based_spp_counts, crown=F, epsilon = 0)
palm_gen_data$r_hat50 <- bd.ms(time=palm_gen_data$age, n=palm_gen_data$tree_based_spp_counts, crown=F, epsilon = 0.5)
palm_gen_data$r_hat90 <- bd.ms(time=palm_gen_data$age, n=palm_gen_data$tree_based_spp_counts, crown=F, epsilon = 0.9)

# intersect with IUCN data
palm_gen_data_old <- palm_gen_data
palm_gen_1 <- palm_gen_data_old[which(palm_gen_data_old$r_hat0 != 0 & palm_gen_data_old$tree_based_spp_counts > 1),]
tmp2_data <- palm_gen_data[which(palm_gen_data$genus %in% names(n_threat)),]
tmp2_data$n_threat<-as.numeric(n_threat[match(tmp2_data$genus,names(n_threat))])
tmp2_data$n_assess<-as.numeric(n_assess[match(tmp2_data$genus,names(n_assess))])
tmp2_data <- tmp2_data[which(tmp2_data$r_hat0 != 0 & tmp2_data$n_assess/tmp2_data$tree_based_spp_counts >= 0.2 & tmp2_data$tree_based_spp_counts > 1),]
tmp2_data <- tmp2_data[match(drop.tip(Palm_GenusTree,setdiff(Palm_GenusTree$tip.label,tmp2_data$genus))$tip.label, tmp2_data$genus),]


sav_vals <- c(# age distribution across tree
              t.test(log(as.numeric(as.character(palm_ages_o_all[[z]]$tip.age))),log(tmp1_data$tip.age), var.equal=T)$p.value,
              t.test(log(as.numeric(as.character(palm_ages_o_all[[z]]$tip.age))),log(tmp1_data$tip.age), var.equal=T)$statistic,
              leveneTest(log(c(tmp1_data$tip.age,as.numeric(as.character(palm_ages_o_all[[z]]$tip.age)))), c(rep("data",dim(tmp1_data)[1]),rep("all",palm_tree_post[[z]]$Nnode+1)))[[3]][1],
              mean(log(tmp1_data$tip.age)),
              sd(log(tmp1_data$tip.age))/sqrt(length(tmp1_data$tip.age)),
              mean(log(as.numeric(as.character(palm_ages_o_all[[z]]$tip.age)))),
              sd(log(as.numeric(as.character(palm_ages_o_all[[z]]$tip.age))))/sqrt(palm_tree_post[[z]]$Nnode+1),

              # age distribution across tree
              t.test(log(palm_gen_1$age),log(tmp2_data$age), var.equal=T)$p.value,
              t.test(log(palm_gen_1$age),log(tmp2_data$age), var.equal=T)$statistic,
              leveneTest(log(c(palm_gen_1$age,tmp2_data$age)), c(rep("data",dim(palm_gen_1)[1]),rep("all",dim(tmp2_data)[1])))[[3]][1],
              mean(log(tmp2_data$age)),
              sd(log(tmp2_data$age))/sqrt(length(tmp2_data$age)),
              mean(log(palm_gen_1$age)),
              sd(log(palm_gen_1$age))/sqrt(length(palm_gen_1$age)),

              # richness distribution across tree
              wilcox.test(log(tmp2_data$tree_based_spp_counts),log(palm_gen_1$tree_based_spp_counts), var.equal=T)$p.value,
              wilcox.test(log(tmp2_data$tree_based_spp_counts),log(palm_gen_1$tree_based_spp_counts), var.equal=T)$statistic[[1]],
              leveneTest(log(c(tmp2_data$tree_based_spp_counts,palm_gen_1$tree_based_spp_counts)), c(rep("data",dim(tmp2_data)[1]),rep("all",dim(palm_gen_1)[1])))[[3]][1],
              mean(log(tmp2_data$tree_based_spp_counts)),
              sd(log(tmp2_data$tree_based_spp_counts))/sqrt(length(tmp2_data$tree_based_spp_counts)),
              mean(log(palm_gen_1$tree_based_spp_counts)),
              sd(log(palm_gen_1$tree_based_spp_counts))/sqrt(length(palm_gen_1$tree_based_spp_counts)),

              # rate distribution across tree
              t.test(log(tmp2_data$r_hat50),log(palm_gen_1$r_hat50), var.equal=T)$p.value,
              t.test(log(tmp2_data$r_hat50),log(palm_gen_1$r_hat50), var.equal=T)$statistic,
              leveneTest(log(c(tmp2_data$r_hat50,palm_gen_1$r_hat50)), c(rep("data",dim(tmp2_data)[1]),rep("all",dim(palm_gen_1)[1])))[[3]][1],
              mean(log(tmp2_data$r_hat50)),
              sd(log(tmp2_data$r_hat50))/sqrt(length(tmp2_data$r_hat50)),
              mean(log(palm_gen_1$r_hat50)),
              sd(log(palm_gen_1$r_hat50))/sqrt(length(palm_gen_1$r_hat50))

             )
return(sav_vals)
})

apply(palm_biases[1+7*0:3,],1,median)
apply(palm_biases[2+7*0:3,],1,median)
apply(palm_biases[4+7*0:3,],1,mean)
apply(palm_biases[5+7*0:3,],1,function(x){sqrt(sum(x^2)/1000)})
apply(palm_biases[6+7*0:3,],1,mean)
apply(palm_biases[7+7*0:3,],1,function(x){sqrt(sum(x^2)/1000)})





################################################################################
###### 15) POWER ANALYSIS FOR THE GENUS-LEVEL RESULTS
################################################################################
p_out_vec <- sapply(seq(25,(length(IUCN_PhytoPhylo_merged$phy$tip)-1),by=25), function(x){
                sapply(1:100, function(y){
                   tmp_tree <- drop.tip(IUCN_PhytoPhylo_merged$phy,setdiff(IUCN_PhytoPhylo_merged$phy$tip.label,IUCN_PhytoPhylo_merged$phy$tip[sample(1:length(IUCN_PhytoPhylo_merged$phy$tip),x,replace=F)]))
                   tmp_tree$node.label <- NULL
                   gls_dat_tmp <- IUCN_ages_sub[match(tmp_tree$tip,IUCN_ages_sub$tip),]
                   rownames(gls_dat_tmp) <- gls_dat_tmp$tip
                   gls_dat_tmp$gls_wts <- 1/sqrt(gls_dat_tmp$n_assess)
                   mod_tmp <- try(gls(I(n_threat/n_assess) ~ log(r_hat50), data= gls_dat_tmp,
                                                                       correlation=corPagel(0.5,tmp_tree,fixed=F),method="ML",
                                                                       weights=varFixed(~1/gls_wts)), silent = T)
                   if(is(mod_tmp,"try-error")) {
                    lam_ints <- seq(0,1,length.out=11)
                    AIC_lam_ints <- sapply(lam_ints, function(z) { AIC(gls(I(n_threat/n_assess) ~ log(r_hat50), data=gls_dat_tmp,correlation=corPagel(z,tmp_tree,fixed=T),method="ML",weights=varFixed(~1/gls_wts))) } )
                    mod_tmp <- gls(I(n_threat/n_assess) ~ log(r_hat50), data=gls_dat_tmp,correlation=corPagel(lam_ints[which.min(AIC_lam_ints)],tmp_tree,fixed=T),method="ML",weights=varFixed(~1/gls_wts))
                   }
                   return(summary(mod_tmp)$tTable['log(r_hat50)',4])
                   })
                 })





################################################################################
###### 16) SENSITIVITY ANALYSIS FOR SAMPLING COVERAGE
################################################################################
# here use entire dataset - don't just subset to those with IUCN classifications - to discover effect of sampling on rate estimates
# find sampling coverage in families
ages1_tmp <- ages1
ages1_tmp$family <- plant_lookup_table[match(ages1_tmp$tip, plant_lookup_table$genus),'family']
ages1_tmp$number.of.accepted.species <- plant_lookup_table[match(ages1_tmp$tip, plant_lookup_table$genus),'number.of.accepted.species']
ages1_tmp$number.of.genera <- as.numeric(genera_in_families[match(ages1_tmp$family, names(genera_in_families))])
ages1_tmp$prop.genera.sampled <- as.numeric(genera_per_family_phylo[match(ages1_tmp$family, names(genera_per_family_phylo))])/ages1_tmp$number.of.genera
ages1_tmp$prop.genera.sampled[is.na(ages1_tmp$prop.genera.sampled)] <- 0

ages1_tmp$r_hat0 <- bd.ms(time=ages1_tmp$tip.age, n=ages1_tmp$number.of.accepted.species, crown=F, epsilon = 0)

ages1_tmp_100cov <- ages1_tmp[which(ages1_tmp$prop.genera.sampled == 1 & ages1_tmp$r_hat0 != 0 & ages1_tmp$number.of.accepted.species > 1),]
PhytoPhylo_ages1_tmp_100cov <- drop.tip(PhytoPhylo_GenusTree,setdiff(PhytoPhylo_GenusTree$tip.label,ages1_tmp_100cov$tip))
PhytoPhylo_ages1_tmp_100cov$node.label <- NULL
ages1_tmp_100cov <- ages1_tmp_100cov[match(PhytoPhylo_ages1_tmp_100cov$tip.label, ages1_tmp_100cov$tip),]

ages1_tmp_100cov_numgen <- table(ages1_tmp_100cov$family)
ages1_tmp_100cov_sub <- ages1_tmp_100cov[ages1_tmp_100cov$family %in% names(ages1_tmp_100cov_numgen[ages1_tmp_100cov_numgen>=4]),]

quantile(with(ages1_tmp_100cov_sub, tapply(number.of.accepted.species, family, sum)),prob=c(0.5,0,1))
quantile(with(ages1_tmp_100cov_sub, tapply(number.of.genera, family, unique)),prob=c(0.5,0,1))


### function to subsample family
ss_fam <- function(fnames, sampling_cov){
df_tmp <- lapply(fnames, function(family_name){
 # first find children in family
 genera_in_family <- plant_lookup_table[which(plant_lookup_table$family == family_name & plant_lookup_table$number.of.accepted.species>0),'genus']
 # confirm they're all in the tree
 if(any(genera_in_family %in% PhytoPhylo_GenusTree$tip == F) == F)
 {
  # find the genera to drop from the tree
  random_genus_samp <- sample(genera_in_family, round((1-sampling_cov)*length(genera_in_family)))
  PhytoPhylo_GenusTree_samp <- drop.tip(PhytoPhylo_GenusTree,random_genus_samp)
  # recalculate ages and rates
  ages1_samp <- calculate_tip_ages(PhytoPhylo_GenusTree_samp)
  ages1_samp[,2] <- as.numeric(as.character(ages1_samp[,2]))
  ages1_samp$family <- plant_lookup_table[match(ages1_samp$tip, plant_lookup_table$genus),'family']
  ages1_samp$number.of.accepted.species <- plant_lookup_table[match(ages1_samp$tip, plant_lookup_table$genus),'number.of.accepted.species']
  ages1_samp$r_hat0 <- bd.ms(time=ages1_samp$tip.age, n=ages1_samp$number.of.accepted.species, crown=F, epsilon = 0)
  merged_samp <- merge(ages1_tmp_100cov_sub[ages1_tmp_100cov_sub$family == family_name,],ages1_samp[ages1_samp$family == family_name,],by='tip')[,c('r_hat0.x','r_hat0.y')]
  return(merged_samp)
   }
 })
 df_tmp <- do.call(rbind,df_tmp)
 #return(cor(df_tmp[,1],df_tmp[,2],method='spearman'))
 return( ((df_tmp[,1] - df_tmp[,2])/df_tmp[,1])*100 )
}

simresults <-  lapply(1:9, function(z){ unlist(sapply(1:50, function(x){ss_fam(fnames = unique(ages1_tmp_100cov_sub$family), sampling_cov = seq(0.9,0.1,by=-0.1)[z])})) })

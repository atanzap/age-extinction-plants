################################################################################
# R code accompanying:
# Tanentzap AJ, Igea J, Johnston MG, Larcombe MJ. 2018. Greater range filling can
# explain why evolutionarily older and slower diversifying plants are less threatened
# by extinction. bioRxiv https://doi.org/10.1101/152215
# See recommendation by Peer Community In Evolutionary Biology at https://dx.doi.org/10.24072/pci.evolbiol.100058
#
# No guarantees provided with this code.
# Prepared by AJ Tanentzap (2 Oct 2018, ajt65 // @ // cam.ac.uk)
#
# Code includes the following sections:
# 1) load packages and relevant functions
# 2) load data for genus-level analysis
# 3) test whether threat status associated with diversification and age for genus-level analysis
# 4) are we taking a bias part of the tol for genus-level analysis?
# 5) load data for species-level analysis
# 6) test association between age and threat status in conifers
# 7) test association between age and range metrics in conifers
# 8) sister species analysis in conifers
# 9) rerun analysis at genus level in conifers
# 10) check bias in the conifer dataset
# 11) test association between age and threat status in palms
# 12) test association between age and range metrics in palms
# 13) sister species analysis in palms
# 14) rerun analysis at genus level in palms
# 15) check bias in the palm dataset

# Data accompanying code available in "./data/" folder of the GitHub repo https://github.com/atanzap/age-extinction-plants

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
  row.names(species.ages)<-tree$tip
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





################################################################################
###### 2) LOAD DATA FOR GENUS-LEVEL ANALYSIS
################################################################################
# read IUCN data
categories.species<-read.csv('./data/export-93437.csv',header=T)
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
IUCN_TPL_Lookup <- lookup_table(as.character(IUCN_TPL_riskunique$Merged),by_species = TRUE,include_counts = TRUE)
rownames(IUCN_TPL_Lookup) <- NULL
IUCN_TPL_Lookup <- IUCN_TPL_Lookup[!duplicated(IUCN_TPL_Lookup$genus),]

# process Qian tree
PhytoPhylo <- read.tree('./data/PhytoPhylo.nex')

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
ii<-sapply(generaPhyndr.2,function(x,y) grep(x,y,fixed=TRUE)[1],y=PhytoPhylo_nooutliers$tip.label)
#drop all but one of every genera
PhytoPhylo_GenusTree<-drop.tip(PhytoPhylo_nooutliers,setdiff(PhytoPhylo_nooutliers$tip.label,PhytoPhylo_nooutliers$tip.label[ii]))
#remove full species name and leave only genera as tip label
PhytoPhylo_GenusTree$tip.label<-sapply(strsplit(PhytoPhylo_GenusTree$tip.label,"_"),function(x) x[1])

# calculate tip ages
ages1 <- calculate_tip_ages(PhytoPhylo_GenusTree)
ages1[,2] <- as.numeric(as.character(ages1[,2]))


# for each genus calculate diversification
IUCN_ages <- merge(ages1,IUCN_TPL_Lookup,by.x="tip",by.y="genus",all.y=T)
IUCN_ages$r_hat0 <- bd.ms(time=IUCN_ages$tip.age, n=IUCN_ages$number.of.accepted.species, crown=F, epsilon = 0)
IUCN_ages$r_hat50 <- bd.ms(time=IUCN_ages$tip.age, n=IUCN_ages$number.of.accepted.species, crown=F, epsilon = 0.5)
IUCN_ages$r_hat90 <- bd.ms(time=IUCN_ages$tip.age, n=IUCN_ages$number.of.accepted.species, crown=F, epsilon = 0.9)


# drop genera that are from families with poor sampling
plant_lookup_table <- plant_lookup(include_counts = TRUE)
genera_in_families <- table(plant_lookup_table[plant_lookup_table$number.of.accepted.species >0,]$family)
IUCN_ages$number.of.genera <- as.numeric(genera_in_families[match(IUCN_ages$family, names(genera_in_families))])
genera_PhytoPhylo <- unique(sapply(strsplit(as.character(PhytoPhylo$tip.label),'_'),function(x) x[1]))
families_PhytoPhylo <- table(plant_lookup_table$family[match(genera_PhytoPhylo,plant_lookup_table$genus)])
IUCN_ages$prop.genera.sampled <- as.numeric(families_PhytoPhylo[match(IUCN_ages$family, names(families_PhytoPhylo))])/IUCN_ages$number.of.genera


# calculate proportion of assessments that are threatened, which accounts for biases in surveying threatened taxa
IUCN_TPL_riskunique_drop <- IUCN_TPL_riskunique[which( is.na(IUCN_TPL_riskunique$Red.List.status) == F & IUCN_TPL_riskunique$Red.List.status != 'DD'),]
threats_per_genus <- with(IUCN_TPL_riskunique_drop, tapply(Red.List.status %in% c("CR","EN","EW","EX","VU"), New.Genus, sum))
assessments_per_genus <- table(IUCN_TPL_riskunique_drop$New.Genus)
IUCN_ages$n_threat <- as.numeric(threats_per_genus[match(IUCN_ages$tip,names(threats_per_genus))])
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
gls_wts <- 1/sqrt(IUCN_ages_sub$n_assess)

# run PGLS weighted by the number of species that are assessed
mod1b <- gls(I(n_threat/n_assess) ~ log(r_hat0), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod2b <- gls(I(n_threat/n_assess) ~ log(r_hat50), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod3b <- gls(I(n_threat/n_assess) ~ log(r_hat90), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod4b <- gls(I(n_threat/n_assess) ~ log(tip.age), data=IUCN_ages_sub, correlation=corPagel(0.5,IUCN_PhytoPhylo_merged$phy,fixed=F),method="ML",weights=varFixed(~1/gls_wts))

cor(fitted(mod1b),with(IUCN_ages_sub,I(n_threat/n_assess)))
cor(fitted(mod2b),with(IUCN_ages_sub,I(n_threat/n_assess)))
cor(fitted(mod3b),with(IUCN_ages_sub,I(n_threat/n_assess)))
cor(fitted(mod4b),with(IUCN_ages_sub,I(n_threat/n_assess)))
with(IUCN_ages_sub,cor(cbind(log(r_hat0),log(r_hat50),log(r_hat90),log(tip.age))))




################################################################################
###### 4) ARE WE TAKING A BIAS PART OF THE TOL FOR GENUS-LEVEL ANALYSIS?
################################################################################
# age distribution across tree
all_ages <- merge(ages1,plant_lookup_table,by.x="tip",by.y="genus",all.y=T)
all_ages$number.of.genera <- as.numeric(genera_in_families[match(all_ages$family, names(genera_in_families))])
all_ages$prop.genera.sampled <- as.numeric(families_PhytoPhylo[match(all_ages$family, names(families_PhytoPhylo))])/all_ages$number.of.genera
all_ages$r_hat0 <- bd.ms(time=all_ages$tip.age, n=all_ages$number.of.accepted.species, crown=F, epsilon = 0)
leveneTest(c(log(IUCN_ages_sub$tip.age),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$tip.age)),
           c(rep("data",dim(IUCN_ages_sub)[1]),rep("all",dim(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$tip.age != 0 & all_ages$number.of.accepted.species > 1),])[1])))
t.test(log(IUCN_ages_sub$tip.age),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$tip.age), var.equal=T)
mean(log(IUCN_ages_sub$tip.age))
sd(log(IUCN_ages_sub$tip.age))/sqrt(length(log(IUCN_ages_sub$tip.age)))
mean(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$tip.age))
sd(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$tip.age))/sqrt(length(log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$tip.age)))

# diversity distribution across tree
leveneTest(c(log(IUCN_ages_sub$number.of.accepted.species),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 50 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species)),
           c(rep("data",dim(IUCN_ages_sub)[1]),rep("all",dim(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$tip.age != 0 & all_ages$number.of.accepted.species > 1),])[1])))
t.test(log(IUCN_ages_sub$number.of.accepted.species),log(all_ages[which(all_ages$prop.genera.sampled >= .6 & all_ages$r_hat0 != 0 & all_ages$number.of.accepted.species > 1),]$number.of.accepted.species), var.equal=T)
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
IUCN_allgens <- as.data.frame(cbind(threats_per_genus, assessments_per_genus))
IUCN_allgens$genus <- as.character(rownames(IUCN_allgens))
rownames(IUCN_allgens) <- NULL
IUCN_allgens <- merge(IUCN_allgens,plant_lookup_table,all.x=T,by='genus')
IUCN_allgens <- IUCN_allgens[which(IUCN_allgens$number.of.accepted.species > 1 & IUCN_allgens$assessments_per_genus/IUCN_allgens$number.of.accepted.species >= 0.2),]
wilcox.test(with(IUCN_allgens,threats_per_genus/assessments_per_genus), with(IUCN_ages_sub, n_threat/n_assess), var.equal=T)
mean(with(IUCN_ages_sub, n_threat/n_assess))
sd(with(IUCN_ages_sub, n_threat/n_assess))/sqrt(length(with(IUCN_ages_sub, n_threat/n_assess)))
mean(with(IUCN_allgens,threats_per_genus/assessments_per_genus))
sd(with(IUCN_allgens,threats_per_genus/assessments_per_genus))/sqrt(length(with(IUCN_allgens,threats_per_genus/assessments_per_genus)))





################################################################################
###### 5) LOAD DATA FOR SPECIES-LEVEL ANALYSIS
################################################################################
#read conifer phylogeny
con_tree <- read.tree('./data/conTree457.tre')

# read the palm phylogeny
palm_tree <- read.tree('./data/palm_tree.tre')

# calculate species ages in conifers, merge with IUCN data and clean up names with TPL
species_ages<-calculate_tip_ages(con_tree)
con_genus<-sapply(as.character(species_ages$tip),function(x) strsplit(x,split='_')[[1]][[1]])
con_species<-sapply(as.character(species_ages$tip),function(x) strsplit(x,split='_')[[1]][[2]])
con_TPL <- TPL(genus=con_genus,species=con_species)
con_TPL$Merged <- with(con_TPL, paste(New.Genus,New.Species,sep='_'))
con_TPL_IUCN <- merge(con_TPL, IUCN_TPL_riskunique_drop[,c('Merged','Red.List.status','Red.List.criteria.version','Year.assessed')], by='Merged', all.x=T)
con_TPL_IUCN$old_species <- with(con_TPL_IUCN, paste(Genus,Species,sep='_'))
species_ages_est <- merge(species_ages, con_TPL_IUCN, by.x = "tip", by.y = "old_species", all.x=T)
species_ages_est$tip <- as.character(species_ages_est$tip)
species_ages_est$tip.age <- as.numeric(as.character(species_ages_est$tip.age))

# calculate species ages in palms and clean up names with TPL
palm_ages_o<-calculate_tip_ages(palm_tree)
palm_genera<-unname(sapply(as.character(palm_ages_o$tip),function(x)strsplit(x,split='_')[[1]][1]))
palm_species<-unname(sapply(as.character(palm_ages_o$tip),function(x)strsplit(x,split='_')[[1]][2]))
palm_TPL<-TPL(genus=palm_genera,species=palm_species)
palm_TPL[,'Merged']<-paste(palm_TPL$New.Genus,palm_TPL$New.Species, sep = '_')
palm_TPL_IUCN <- merge(palm_TPL, IUCN_TPL_riskunique_drop[,c('Merged','Red.List.status','Red.List.criteria.version','Year.assessed')], by='Merged', all.x=T)
palm_TPL_IUCN$old_species <- with(palm_TPL_IUCN, paste(Genus,Species,sep='_'))
palm_ages <- merge(palm_ages_o, palm_TPL_IUCN, by.x = "tip", by.y = "old_species", all.x=T)
palm_ages$tip <- as.character(palm_ages$tip)
palm_ages$tip.age <- as.numeric(as.character(palm_ages$tip.age))

#remove synonyms
species_ages_est_a <- species_ages_est[species_ages_est$Taxonomic.status=='Accepted',]
palm_ages <- palm_ages[palm_ages$Taxonomic.status=='Accepted',]
palm_dups <- palm_ages[duplicated(palm_ages$Merged),]$Merged
palm_ages <- palm_ages[-which(palm_ages$Merged %in% palm_dups & (with(palm_ages,paste(Genus,Species,sep="_")) %in% palm_ages$Merged == F)),]

# simplify risk status
levels(species_ages_est_a$Red.List.status)[levels(species_ages_est_a$Red.List.status)%in%c("CR","EN","VU")] <- "Threatened"
levels(species_ages_est_a$Red.List.status)[levels(species_ages_est_a$Red.List.status)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
species_ages_est_a <- species_ages_est_a[species_ages_est_a$Red.List.status %in% c('EW','EX','DD') == F,]
species_ages_est_a$Red.List.status <- as.numeric(species_ages_est_a$Red.List.status=="Threatened")
levels(palm_ages$Red.List.status)[levels(palm_ages$Red.List.status)%in%c("CR","EN","VU")] <- "Threatened"
levels(palm_ages$Red.List.status)[levels(palm_ages$Red.List.status)%in%c("LC","LR/cd","LR/lc","LR/nt","NT")] <- "Lower risk"
palm_ages <- palm_ages[palm_ages$Red.List.status %in% c('EW','EX','DD') == F,]
palm_ages$Red.List.status <- as.numeric(palm_ages$Red.List.status=="Threatened")





################################################################################
###### 6) TEST ASSOCIATION BETWEEN AGE AND THREAT STATUS IN CONIFERS
################################################################################
IUCN_Qian_data_merged <- comparative.data(phy=con_tree,data=species_ages_est_a[,c('tip','tip.age','Red.List.status')],names.col=tip,vcv=T)
rownames(IUCN_Qian_data_merged$data) <- IUCN_Qian_data_merged$phy$tip
m4a <- phyloglm(Red.List.status~scale(log(tip.age)),phy=IUCN_Qian_data_merged$phy,data=IUCN_Qian_data_merged$data,boot=100,method="logistic_MPLE", log.alpha.bound = 10)




################################################################################
###### 7) TEST ASSOCIATION BETWEEN AGE AND RANGE METRICS IN CONIFERS
################################################################################
# load distribution data
conifer_ranges <- read.csv('./data/conifer_rangeSize_table.csv', header = T)
species_data_sub <- merge(species_ages_est_a,conifer_ranges,by.x='tip',by.y='species')
species_data_sub$gbif_recs <- species_data_sub$d.tpos + species_data_sub$d.fneg
species_data_sub$log_rfilled_g <- with(species_data_sub, qlogis(EA_0.25_cells/global.count))

con_data_merged_all <- comparative.data(phy=con_tree,data=species_data_sub[,c('tip','tip.age','Red.List.status','log_rfilled_g','gbif_recs','EA_0.25_cells','global.count')],names.col=tip,vcv=T)
m5a <- pgls(log(global.count)~scale(log(tip.age))*as.factor(Red.List.status), con_data_merged_all, lambda='ML')
m5b <- pgls(log(gbif_recs)~log(tip.age)*as.factor(Red.List.status), con_data_merged_all, lambda='ML')





################################################################################
###### 8) SISTER SPECIES ANALYSIS IN CONIFERS
################################################################################
# find sisters in full tree and then use those with contrasting extinction status
sister_extract <- cbind(con_tree$tip,as.character(unlist(sapply(con_tree$tip,getSisters,tree=con_tree,mode='label'))))
sister_extract <- sister_extract[!duplicated(t(apply(sister_extract, 1, sort))), ]
sister_extract <- sister_extract[apply(sister_extract,1,function(x){sum(is.na(as.numeric(x)))})==2,]
sister_extract <- cbind(sister_extract, as.data.frame(t(apply(sister_extract, 1, function(x){  con_data_merged_all$data[match(x, rownames(con_data_merged_all$data)),]$Red.List.status  }))))
colnames(sister_extract) <- c('Sp1','Sp2','Stat1','Stat2')
sister_extract <- sister_extract[which(is.na(sister_extract$Stat1) == F & is.na(sister_extract$Stat2) == F),]
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
tru_cors <- sapply(which(colnames(sister_extract_dat)=='diff_PotRange'):ncol(sister_extract_dat), function(z){ cor(sister_extract_dat[,z],log(sister_extract_dat$Time)) })
names(tru_cors) <- 'diff_PotRange'

sapply(1:ncol(ran_cors),function(x){
 hist(ran_cors[,x],breaks=20, main = paste(colnames(sister_extract_dat)[which(colnames(sister_extract_dat)=='diff_PotRange'):ncol(sister_extract_dat)][x], " = ", sum(tru_cors[x] > ran_cors[,x])/1000))
 abline(v=tru_cors[x],lwd=2,col='red')
 })





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
plant_lookup_table <- plant_lookup(include_counts = TRUE)
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
con_TPL_IUCN2$Red.List.status <- as.numeric(con_TPL_IUCN2$Red.List.status=="Threatened")
n_threat <- with(con_TPL_IUCN2, tapply(Red.List.status,Genus[drop=T],sum))
n_assess <- with(con_TPL_IUCN2, tapply(Red.List.status,Genus[drop=T],length))

con_gen_data <- con_gen_data[which(con_gen_data$genus %in% names(n_threat)),]
con_gen_data$n_threat<-as.numeric(n_threat[match(con_gen_data$genus,names(n_threat))])
con_gen_data$n_assess<-as.numeric(n_assess[match(con_gen_data$genus,names(n_assess))])

# prepare data for PGLS
con_gen_data_sub <- con_gen_data[which(con_gen_data$r_hat0 != 0 & con_gen_data$n_assess/con_gen_data$number.of.accepted.species >= 0.2 & con_gen_data$number.of.accepted.species > 1),]
Con_GenusTree_IUCNsub <- drop.tip(Con_GenusTree,setdiff(Con_GenusTree$tip.label,con_gen_data_sub$genus))
con_gen_data_sub <- con_gen_data_sub[match(Con_GenusTree_IUCNsub$tip.label, con_gen_data_sub$genus),]

# re-run PGLS weighted by the number of species that are assessed
gls_wts <- 1/sqrt(con_gen_data_sub$n_assess)
mod1b <- gls(I(n_threat/n_assess) ~ log(r_hat0),  data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod2b <- gls(I(n_threat/n_assess) ~ log(r_hat50), data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod3b <- gls(I(n_threat/n_assess) ~ log(r_hat90), data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
mod4b <- gls(I(n_threat/n_assess) ~ log(age), data=con_gen_data_sub, correlation=corPagel(0.5,Con_GenusTree_IUCNsub,fixed=F),method="ML",weights=varFixed(~1/gls_wts))
round(cor(cbind(fitted(mod1b),fitted(mod2b),fitted(mod3b),fitted(mod4b),with(con_gen_data_sub,I(n_threat/n_assess)))),2)





################################################################################
###### 10) CHECK BIAS IN THE CONIFER DATASET
################################################################################
# species age distribution across tree
con_TPL2 <- con_TPL
con_TPL2$old_species <- with(con_TPL2, paste(Genus,Species,sep='_'))
species_ages_allcons <- merge(species_ages, con_TPL2, by.x = "tip", by.y = "old_species", all.x=T)
tmp_con <- species_ages_allcons[species_ages_allcons$Taxonomic.status=='Accepted',]
leveneTest(log(c(IUCN_Qian_data_merged$data$tip.age,as.numeric(as.character(tmp_con$tip.age)))), c(rep("data",dim(IUCN_Qian_data_merged$data)[1]),rep("all",dim(tmp_con)[1])))
t.test(log(IUCN_Qian_data_merged$data$tip.age),log(as.numeric(as.character(tmp_con$tip.age))), var.equal=T)
mean(log(IUCN_Qian_data_merged$data$tip.age))
sd(log(IUCN_Qian_data_merged$data$tip.age))/sqrt(length(IUCN_Qian_data_merged$data$tip.age))
mean(log(as.numeric(as.character(tmp_con$tip.age))))
sd(log(as.numeric(as.character(tmp_con$tip.age))))/sqrt(length(tmp_con$tip.age))





################################################################################
###### 11) TEST ASSOCIATION BETWEEN AGE AND THREAT STATUS IN PALMS
################################################################################
palms_IUCN_drop <- palm_ages[is.na(palm_ages$Red.List.status)==F,]
palms_dropped_tree <- drop.tip(palm_tree, setdiff(palm_tree$tip, palms_IUCN_drop$tip) )
palms_dropped_tree$tip.label <- palms_IUCN_drop[match(palms_dropped_tree$tip.label, palms_IUCN_drop$tip),'Merged']
palms_IUCN_drop <- palms_IUCN_drop[match(palms_dropped_tree$tip.label,palms_IUCN_drop$Merged),]
rownames(palms_IUCN_drop) <- palms_IUCN_drop$Merged
mp1 <- phyloglm(Red.List.status~scale(log(tip.age)),phy=palms_dropped_tree,data=palms_IUCN_drop,method="logistic_MPLE", log.alpha.bound = 10)





################################################################################
###### 12) TEST ASSOCIATION BETWEEN AGE AND RANGE METRICS IN PALMS
################################################################################
# load TTR and GBIF raw data
palm_stats <- read.csv('./data/palm_stats.csv')
palm_ages_all <- merge(palm_ages,palm_stats,by.x='Merged',by.y='sp')
palm_ages_all$gbif_recs <- palm_ages_all$d.tpos + palm_ages_all$d.fneg
palm_ages_all$log_rfilled_g <- with(palm_ages_all, qlogis(EA_0.25_cells/global.count))
palms_IUCN_drop_range <- merge(palms_IUCN_drop,palm_ages_all[,c('Merged','EA_0.25_cells','global.count','gbif_recs','log_rfilled_g')],by='Merged')
palm_data_merged <- comparative.data(phy=palm_tree,data=palms_IUCN_drop_range[,c('tip','tip.age','Red.List.status','log_rfilled_g','gbif_recs','EA_0.25_cells','global.count')],names.col=tip,vcv=T)

mp2a <- pgls(log(global.count)~log(tip.age)*as.factor(Red.List.status), palm_data_merged, lambda='ML')
mp2b <- pgls(log(gbif_recs)~log(tip.age)*as.factor(Red.List.status), palm_data_merged, lambda='ML')





################################################################################
###### 13) SISTER SPECIES ANALYSIS IN PALMS
################################################################################
# find sisters with contrasting extinction status
sister_extract_p <- cbind(palm_data_merged$phy$tip,as.character(unlist(sapply(palm_data_merged$phy$tip,getSisters,tree=palm_tree,mode='label'))))
sister_extract_p <- sister_extract_p[!duplicated(t(apply(sister_extract_p, 1, sort))), ]
sister_extract_p <- sister_extract_p[apply(sister_extract_p,1,function(x){sum(is.na(as.numeric(x)))})==2,]
sister_extract_p <- cbind(sister_extract_p, as.data.frame(t(apply(sister_extract_p, 1, function(x){  palms_IUCN_drop_range[match(x, with(palms_IUCN_drop_range,paste(Genus,Species,sep="_"))),]$Red.List.status  }))))
colnames(sister_extract_p) <- c('Sp1','Sp2','Stat1','Stat2')
sister_extract_p <- sister_extract_p[which(is.na(sister_extract_p$Stat1) == F & is.na(sister_extract_p$Stat2) == F),]
sister_extract_p_nomatch <- sister_extract_p[which(sister_extract_p$Stat1 != sister_extract_p$Stat2),]

# extract characteristics for data frame of sisters
rownames(palms_IUCN_drop_range) <- palms_IUCN_drop_range$tip
sister_extract_p_nomatch[,1] <- as.character(sister_extract_p_nomatch[,1])
sister_extract_p_nomatch[,2] <- as.character(sister_extract_p_nomatch[,2])
sister_extract_dat_p <- as.data.frame(t( sapply(1:nrow(sister_extract_p_nomatch), function(x){sis_chars(sister_extract_p_nomatch[x,], dataf = palms_IUCN_drop_range, sampling = F)}) ))
colnames(sister_extract_dat_p) <- c('Sp1','Sp2','Time','diff_PotRange')
sister_extract_dat_p[,3:ncol(sister_extract_dat_p)] <- sapply(sister_extract_dat_p[,3:ncol(sister_extract_dat_p)], function(x){ as.numeric(as.character(x))})

# find sisters that are both threatened and sample from them - we can't do non-threatened because there are only 2 options
sister_extract_match_p <- sister_extract_p[which(sister_extract_p$Stat1 == sister_extract_p$Stat2 & sister_extract_p$Stat2 == 1),]
sister_extract_match_p[,1] <- as.character(sister_extract_match_p[,1])
sister_extract_match_p[,2] <- as.character(sister_extract_match_p[,2])

ran_cors_p <- lapply(1:1000, function(y){
   newdata <- sister_extract_match_p[sample(1:nrow(sister_extract_match_p), nrow(sister_extract_p_nomatch), replace=F), ]
   sister_extract_dat2 <- as.data.frame(t( sapply(1:nrow(sister_extract_p_nomatch), function(x){sis_chars(newdata[x,], dataf = palms_IUCN_drop_range, sampling = T)}) ))
   colnames(sister_extract_dat2) <- c('Sp1','Sp2','Time','diff_PotRange')
   sister_extract_dat2[,3:ncol(sister_extract_dat2)] <- sapply(sister_extract_dat2[,3:ncol(sister_extract_dat2)], function(x){ as.numeric(as.character(x))})
   output <- sapply(which(colnames(sister_extract_dat2)=='diff_PotRange'):ncol(sister_extract_dat2), function(z){ cor(sister_extract_dat2[,z],log(sister_extract_dat2$Time)) })
   return(output)
  })
ran_cors_p <- do.call("rbind",ran_cors_p)
tru_cors_p <- sapply(which(colnames(sister_extract_dat_p)=='diff_PotRange'):ncol(sister_extract_dat_p), function(z){ cor(sister_extract_dat_p[,z],log(sister_extract_dat_p$Time)) })
names(tru_cors_p) <- c('diff_PotRange')

sapply(1:ncol(ran_cors_p),function(x){
 hist(ran_cors_p[,x],breaks=20,main=paste(colnames(sister_extract_dat_p)[which(colnames(sister_extract_dat_p)=='diff_PotRange'):ncol(sister_extract_dat_p)][x]," = ",sum(tru_cors_p[x] > ran_cors_p[,x])/1000))
 abline(v=tru_cors_p[x],lwd=2,col='red')
})





################################################################################
###### 14) RERUN ANALYSIS AT GENUS LEVEL IN PALMS
################################################################################
# filter tree
generaPhyndr <- unique(sapply(strsplit(as.character(palms_dropped_tree$tip.label),'_'),function(x) x[1]))
generaPhyndr.2 <- sapply(generaPhyndr, function(x) x<- paste(x,'_',sep=''))
ii<-sapply(generaPhyndr.2,function(x,y) grep(x,y,fixed=TRUE)[1],y=palms_dropped_tree$tip.label)
#drop all but one of every genera
Palm_GenusTree<-drop.tip(palms_dropped_tree,setdiff(palms_dropped_tree$tip.label,palms_dropped_tree$tip.label[ii]))
#remove full species name and leave only genera as tip label
Palm_GenusTree$tip.label<-sapply(strsplit(Palm_GenusTree$tip.label,"_"),function(x) x[1])

#assemble dataframe
palm_gen_data <- as.data.frame(cbind( names(with(palm_ages_o, tapply(as.numeric(as.character(tip.age)), sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max))),
                                      as.numeric(with(palm_ages_o, tapply(as.numeric(as.character(tip.age)), sapply(strsplit(as.character(tip),'_'),function(x) x[1]), max)))
                    ))
palm_gen_data[,1] <- as.character(palm_gen_data[,1])
palm_gen_data[,2] <- as.numeric(as.character(palm_gen_data[,2]))
colnames(palm_gen_data) <- c("genus","age")
plant_lookup_table <- plant_lookup(include_counts = TRUE)
palm_gen_data <- merge(palm_gen_data,plant_lookup_table,all.x=T,by="genus")
n_threat <- with(palms_IUCN_drop, tapply(Red.List.status,Genus[drop=T],sum))
n_assess <- with(palms_IUCN_drop, tapply(Red.List.status,Genus[drop=T],length))

# derive spp counts from tree rather than TPL (some changes don't appear in TPL)
tree_based_spp_counts <- table(sapply(strsplit(as.character(palm_tree$tip.label),'_'),function(x) x[1]))
palm_gen_data$tree_based_spp_counts <- tree_based_spp_counts[match(names(tree_based_spp_counts),palm_gen_data$genus)]
#with(palm_gen_data, plot(number.of.accepted.species,tree_based_spp_counts))

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
gls_wts <- 1/sqrt(palm_gen_data_sub$n_assess)

mod1p <- gls(I(n_threat/n_assess) ~ log(r_hat0),  data=palm_gen_data_sub, correlation=corPagel(0,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))
mod2p <- gls(I(n_threat/n_assess) ~ log(r_hat50), data=palm_gen_data_sub, correlation=corPagel(0,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))
mod3p <- gls(I(n_threat/n_assess) ~ log(r_hat90), data=palm_gen_data_sub, correlation=corPagel(0,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))
mod4p <- gls(I(n_threat/n_assess) ~ log(age),     data=palm_gen_data_sub, correlation=corPagel(1,Palm_GenusTree_IUCNsub,fixed=T),method="ML",weights=varFixed(~1/gls_wts))
round(cor(cbind(fitted(mod1p),fitted(mod2p),fitted(mod3p),fitted(mod4p),with(palm_gen_data_sub,I(n_threat/n_assess)))),2)





################################################################################
###### 15) CHECK BIAS IN THE PALM DATASET
################################################################################
# age distribution across tree
leveneTest(log(c(palms_IUCN_drop$tip.age,as.numeric(as.character(palm_ages_o$tip.age)))), c(rep("data",dim(palms_IUCN_drop)[1]),rep("all",palm_tree$Nnode+1)))
t.test(log(as.numeric(as.character(palm_ages_o$tip.age))),log(palms_IUCN_drop$tip.age), var.equal=T)
mean(log(palms_IUCN_drop$tip.age))
sd(log(palms_IUCN_drop$tip.age))/sqrt(length(palms_IUCN_drop$tip.age))
mean(log(as.numeric(as.character(palm_ages_o$tip.age))))
sd(log(as.numeric(as.character(palm_ages_o$tip.age))))/sqrt(palm_tree$Nnode+1)

# age distribution across tree
palm_tmp1 <- palm_gen_data_old[which(palm_gen_data_old$r_hat0 != 0 & palm_gen_data_old$tree_based_spp_counts > 1),]
leveneTest(log(c(palm_tmp1$age,palm_gen_data_sub$age)), c(rep("data",dim(palm_tmp1)[1]),rep("all",dim(palm_gen_data_sub)[1])))
t.test(log(palm_tmp1$age),log(palm_gen_data_sub$age), var.equal=T)
mean(log(palm_tmp1$age))
sd(log(palm_tmp1$age))/sqrt(length(palm_tmp1$age))
mean(log(palm_gen_data_sub$age))
sd(log(palm_gen_data_sub$age))/sqrt(length(palm_gen_data_sub$age))

# richness distribution across tree
leveneTest(log(c(palm_tmp1$tree_based_spp_counts,palm_gen_data_sub$tree_based_spp_counts)), c(rep("data",dim(palm_tmp1)[1]),rep("all",dim(palm_gen_data_sub)[1])))
wilcox.test(log(palm_tmp1$tree_based_spp_counts),log(palm_gen_data_sub$tree_based_spp_counts), var.equal=T)
mean(log(palm_tmp1$tree_based_spp_counts))
sd(log(palm_tmp1$tree_based_spp_counts))/sqrt(length(palm_tmp1$tree_based_spp_counts))
mean(log(palm_gen_data_sub$tree_based_spp_counts))
sd(log(palm_gen_data_sub$tree_based_spp_counts))/sqrt(length(palm_gen_data_sub$tree_based_spp_counts))

# rate distribution across tree
leveneTest(log(c(palm_tmp1$r_hat50,palm_gen_data_sub$r_hat50)), c(rep("data",dim(palm_tmp1)[1]),rep("all",dim(palm_gen_data_sub)[1])))
t.test(log(palm_tmp1$r_hat50),log(palm_gen_data_sub$r_hat50), var.equal=T)
mean(log(palm_tmp1$r_hat50))
sd(log(palm_tmp1$r_hat50))/sqrt(length(palm_tmp1$r_hat50))
mean(log(palm_gen_data_sub$r_hat50))
sd(log(palm_gen_data_sub$r_hat50))/sqrt(length(palm_gen_data_sub$r_hat50))

##
## R source code to accompany <publication>, last updated 14 August 2018.
## Please contact Ava Hoffman (avamariehoffman@gmail.com) with questions.
##
## Ensure you have dependent files:
## 
## "SpComp\ -\ Genetic\ Diversity\ -\ KNZ\ 2014\ COMPLETED.csv"
## "SpComp\ -\ Genetic\ Diversity\ -\ SGS\ 2014\ COMPLETED.csv"
##
## If you found this code useful, please use the citation below:
##
#############################################################################################
## Set working directory

wd <-
  "/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/phylogenetic-diversity-pilot"
setwd(wd)

KNZ <-  read.csv("SpComp\ -\ Genetic\ Diversity\ -\ KNZ\ 2014\ COMPLETED.csv", header = T)
SGS <- read.csv("SpComp\ -\ Genetic\ Diversity\ -\ SGS\ 2014\ COMPLETED.csv", header = T)
ILL <- read.csv("spComp_GeneticDiversity_illinois_compiled_blk1removed.csv", header=T)


###########################################################################################
## Load packages
library(ape)
library(tidyverse)
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")
library(ggtree)
library(ggplot2)
library(reshape)
library(picante)
library(diverse)
library(codyn)

###########################################################################################
###########################################################################################
##
## SGS 
##
###########################################################################################
###########################################################################################

## need to modify input file to look nexus file
## Add this to the beginning:
## #NEXUS
## 
## begin trees;
## trees tree1 = 
## And this to the end:
## end;

community_SGS <- read.nexus("/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/phylogenetic-diversity-pilot/SGS/Cipres_Data/RAxML_bestTree.result")
community_SGS$edge.length <- round(community_SGS$edge.length, digits = 3)
community_SGS.drop <- (subtrees(community_SGS)[[2]]) # don't want the outgroup
## couldn't figure out a better way to detect base tree other than running subtrees(community_SGS) and browsing...
ggtree(community_SGS.drop) +
  geom_tiplab(size=7) +
  #geom_text(aes(label=node), hjust=-1, vjust=0.1) +
  geom_text(aes(label=branch.length), hjust=1, vjust=0.1, size=4) +
  geom_hilight(node=35, fill="darkgreen", alpha=.3) +
  geom_hilight(node=43, fill="steelblue", alpha=.3) +
  xlim(NA, 0.5)
ggsave(width=30,height=15, filename="SGS/Distances_toscale.jpg")

###########################################################################################
## process species comp data
SGS_df <- SGS[,1:5] # just keep the first five columns
# remove unknowns
SGS_df <- SGS_df[!(SGS_df$species == "Unk_6-_smallaster"),]
SGS_df <- SGS_df[!(SGS_df$species == "bare_ground"),]
SGS_df <- SGS_df[!(SGS_df$species == "cryptogramic"),]
SGS_df$species <- as.character(SGS_df$species)
# a few have antiquated names
SGS_df[SGS_df == "Buchloe_dactyloides"] <- "Bouteloua_dactyloides"
SGS_df[SGS_df == "Psoralea_tenuiflora"] <- "Psoralidium_tenuiflorum"
    ## how many unique species?
    print("number of species");length(sort(unique(SGS_df$species))) ## 33 spp
    sort(unique(SGS_df$species))
clean_SGS <- data.frame()
for(i in 1:3){
  for(j in 1:10){
    PlotID <- SGS_df[(SGS_df$plot == j),]
    PlotID <- PlotID[(PlotID$Block == i),]
    PlotID <- PlotID[order(PlotID[,'species'],-PlotID[,'X.cover']),] ## sort so largest percent cover is on top
    PlotID <- PlotID[!duplicated(PlotID$species),] ## remove any species duplicated with lower percent cover
    print(PlotID)
    clean_SGS <- rbind(clean_SGS,PlotID)
  }
}
write.csv(clean_SGS, "Clean_spp_comp_SGS.csv")
# remove date
clean_SGS_wide <- clean_SGS[,-1]
# reshape to wide format
clean_SGS_wide <- reshape(clean_SGS_wide, direction = "wide", timevar = "species", idvar = c("Block","plot"))
# fill NAs with zeros
clean_SGS_wide[is.na(clean_SGS_wide)] <- 0
# clean up species names
names(clean_SGS_wide) <- gsub("X.cover.", "", names(clean_SGS_wide))
# remove block and plot column
clean_SGS_wide <- clean_SGS_wide[,-c(1:2)]
# correct rownames
rownames(clean_SGS_wide) <- seq(1,30,1)

## calculate Faith's PD (see picante package)
tree_SGS <- prune.sample(clean_SGS_wide,community_SGS)
phylo.div_SGS <- pd(clean_SGS_wide,community_SGS, include.root = F)

## Shannons + Simpsons's
div_SGS_Shan <- community_diversity(clean_SGS, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "Shannon")
div_SGS_Simp <- community_diversity(clean_SGS, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "InverseSimpson")

## Berger Parker
bp.index_SGS <- diversity(as.matrix(clean_SGS_wide), type = "berger-parker")
# reorder
bp.index_SGS <- bp.index_SGS[order(as.numeric(row.names(bp.index_SGS))),]

## evenness  - use codyn package
evenness_SGS_EQ <- community_structure(clean_SGS, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "EQ")
evenness_SGS_Evar <- community_structure(clean_SGS, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "Evar")

# all metrics
phylo.div_SGS 
div_SGS_Shan
div_SGS_Simp
bp.index_SGS
evenness_SGS_EQ
evenness_SGS_Evar

## concatenate
total.div_SGS <- merge(div_SGS_Shan,div_SGS_Simp)
total.div_SGS <- merge(total.div_SGS,evenness_SGS_EQ)
total.div_SGS <- merge(total.div_SGS,evenness_SGS_Evar)
total.div_SGS <- total.div_SGS[order(total.div_SGS$Block,total.div_SGS$plot),]
total.div_SGS <- cbind(total.div_SGS, bp.index_SGS, phylo.div_SGS)
rownames(total.div_SGS) <- seq(1,30,1)
# sanity check = are richness and SR the same
total.div_SGS[,c("richness","SR")]
# looks good
write.csv(total.div_SGS, "DIVERSITY_sgs.csv")

###########################################################################################
## import biomass data
biomassdata = read.csv("Biomass_data_all.csv",header=TRUE)
biomassdata_SGS = biomassdata[(biomassdata$Site == "SGS"),]
phylodata_SGS = read.csv("DIVERSITY_sgs.csv",header=TRUE)
phylodata_SGS <- phylodata_SGS[,-1]

nrow(biomassdata_SGS)
nrow(phylodata_SGS)

completedata_SGS <-  merge(biomassdata_SGS, phylodata_SGS)
# remove irrelevant columns
completedata_SGS <-  subset(completedata_SGS, select=-c(ANGE,Woody))
# rename column
completedata_SGS$Graminoids <-  completedata_SGS$Grass 
# total cactus is the total without cactus mass
completedata_SGS$Total <-  completedata_SGS$total.cactus ## exclude cactus mass..
#completedata_SGS$pd.obs.per.spp <- completedata_SGS$pd.obs / completedata_SGS$ntaxa ## calculate phylogenetic distance per species
#just want a subset of variables for plotting
df.part_SGS <- subset(completedata_SGS, select=c(Block,plot,BOGR.BOHI,Forbs,Graminoids,Total,PD,richness,berger.parker.D,Evar,Shannon,InverseSimpson))
mcd_SGS <- melt(df.part_SGS, id.vars = c("Block","plot","PD","richness","berger.parker.D","Evar","Shannon","InverseSimpson"))
mcd_SGS <- melt(mcd_SGS, id.vars = c("Block","plot","variable","value"))
colnames(mcd_SGS) <- c("Block","plot","bio_type","biomass","metric_type","metric")

#calculate p.values
p.df <- data.frame(0,0,0)
for(i in 1:4){
  for(j in 1:6){
    bt=levels(mcd_SGS$bio_type)[i]
    mt=levels(mcd_SGS$metric_type)[j]
    p.val=coef(summary(lm(biomass ~ metric, data = mcd_SGS[ mcd_SGS$bio_type==bt & mcd_SGS$metric_type==mt, ])))[2,4]
    dats <- c(bt,mt,p.val)
    p.df <- rbind(p.df,dats)
  }
}
colnames(p.df) <- c("bio_type","metric_type","p.val")
p.df <- p.df[-1,]

mcd_SGS_pval <- merge(mcd_SGS,p.df)

## plot all diversity metrics and biomass
ggplot(data = mcd_SGS_pval,
       aes(x=metric,y=biomass)) +
  facet_wrap(metric_type~bio_type, scales = "free", nrow = 7 , ncol = 4)+
  geom_smooth(data=subset(mcd_SGS_pval,  p.val < 0.05), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, aes(color=as.numeric(p.val))) +
  geom_point() +
  theme_classic() +
  xlab("Diversity metric") +
  ylab("Biomass") +
  guides(color=guide_legend(title="p-value"))
ggsave(filename="Figures/Biomass_v_Diversity_SGS.jpg",height=14,width=10)

###########################################################################################
###########################################################################################
##
## Konza
##
###########################################################################################
###########################################################################################
## phylogenetic tree

community_KNZ <- read.nexus("/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/phylogenetic-diversity-pilot/Konza/Cipres_Data/RAxML_bestTree.result")
community_KNZ$edge.length <- round(community_KNZ$edge.length, digits = 3)
community_KNZ.drop <- (subtrees(community_KNZ)[[2]]) # don't want the outgroup
ggtree(community_KNZ.drop) +
  geom_tiplab(size=7) +
  #geom_text(aes(label=node), hjust=-1, vjust=0.1) +
  geom_text(aes(label=branch.length), hjust=1, vjust=0.1, size=4) +
  geom_hilight(node=51, fill="darkgreen", alpha=.3) +
  geom_hilight(node=65, fill="steelblue", alpha=.3) +
  xlim(NA, 0.68)
ggsave(width=30,height=18, filename="Konza/Distances_toscale.jpg")

###########################################################################################
## process species comp data

KNZ_df <- KNZ
## how many unique species?
print("number of species");length(sort(unique(KNZ_df$species))) ## 49 spp
sort(unique(KNZ_df$species))
clean_KNZ <- data.frame()
for(i in 1:3){
  for(j in 1:10){
    PlotID <- KNZ_df[(KNZ_df$plot == j),]
    PlotID <- PlotID[(PlotID$Block == i),]
    PlotID <- PlotID[order(PlotID[,'species'],-PlotID[,'X.cover']),] ## sort so largest percent cover is on top
    PlotID <- PlotID[!duplicated(PlotID$species),] ## remove any species duplicated with lower percent cover
    print(PlotID)
    clean_KNZ <- rbind(clean_KNZ,PlotID)
  }
}
write.csv(clean_KNZ, "Clean_spp_comp_KNZ.csv")
clean_KNZ_wide <- clean_KNZ[,-1]
clean_KNZ_wide <- reshape(clean_KNZ_wide, direction = "wide", timevar = "species", idvar = c("Block","plot"))
clean_KNZ_wide[is.na(clean_KNZ_wide)] <- 0
names(clean_KNZ_wide) <- gsub("X.cover.", "", names(clean_KNZ_wide))
clean_KNZ_wide <- clean_KNZ_wide[,-c(1:2)]
rownames(clean_KNZ_wide) <- seq(1,30,1)

## calculate Faith's PD (see picante package)
tree_KNZ <- prune.sample(clean_KNZ_wide,community_KNZ)
phylo.div_KNZ <- pd(clean_KNZ_wide,community_KNZ, include.root = F)

## Shannons + Simpsons's
div_KNZ_Shan <- community_diversity(clean_KNZ, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "Shannon")
div_KNZ_Simp <- community_diversity(clean_KNZ, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "InverseSimpson")

## Berger Parker
bp.index_KNZ <- diversity(as.matrix(clean_KNZ_wide), type = "berger-parker")
bp.index_KNZ <- bp.index_KNZ[order(as.numeric(row.names(bp.index_KNZ))),]

## evenness
evenness_KNZ_EQ <- community_structure(clean_KNZ, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "EQ")
evenness_KNZ_Evar <- community_structure(clean_KNZ, replicate.var = c("plot","Block"), abundance.var = "X.cover", metric = "Evar")

## concatenate
total.div_KNZ <- merge(div_KNZ_Shan,div_KNZ_Simp)
total.div_KNZ <- merge(total.div_KNZ,evenness_KNZ_EQ)
total.div_KNZ <- merge(total.div_KNZ,evenness_KNZ_Evar)
total.div_KNZ <- total.div_KNZ[order(total.div_KNZ$Block,total.div_KNZ$plot),]
total.div_KNZ <- cbind(total.div_KNZ, bp.index_KNZ, phylo.div_KNZ)
rownames(total.div_KNZ) <- seq(1,30,1)
# sanity check = are richness and SR the same
total.div_KNZ[,c("richness","SR")]
# 
write.csv(total.div_KNZ, "DIVERSITY_KNZ.csv")

###########################################################################################
## import biomass data
biomassdata = read.csv("Biomass_data_all.csv",header=TRUE)
biomassdata_KNZ = biomassdata[(biomassdata$Site == "KNZ"),]
phylodata_KNZ = read.csv("DIVERSITY_KNZ.csv",header=TRUE)
phylodata_KNZ <- phylodata_KNZ[,-1]

nrow(biomassdata_KNZ)
nrow(phylodata_KNZ)

completedata_KNZ <-  merge(biomassdata_KNZ, phylodata_KNZ)
# remove irrelevant columns
completedata_KNZ <-  subset(completedata_KNZ, select=-c(BOGR.BOHI,Cactus,total.cactus))
completedata_KNZ$Graminoids <-  completedata_KNZ$Grass + completedata_KNZ$ANGE
completedata_KNZ$Total <- completedata_KNZ$total
#completedata_KNZ$pd.obs.per.spp <- completedata_KNZ$pd.obs / completedata_KNZ$ntaxa ## calculate phylogenetic distance per species
df.part_KNZ <- subset(completedata_KNZ, select=c(Block,plot,ANGE,Forbs,Graminoids,Woody,Total,PD,richness,berger.parker.D,Evar,Shannon,InverseSimpson))
mcd_KNZ <- melt(df.part_KNZ, id.vars = c("Block","plot","PD","richness","berger.parker.D","Evar","Shannon","InverseSimpson"))
mcd_KNZ <- melt(mcd_KNZ, id.vars = c("Block","plot","variable","value"))
colnames(mcd_KNZ) <- c("Block","plot","bio_type","biomass","metric_type","metric")

#calculate p.values
p.df <- data.frame(0,0,0)
for(i in 1:5){
  for(j in 1:6){
    bt=levels(mcd_KNZ$bio_type)[i]
    mt=levels(mcd_KNZ$metric_type)[j]
    p.val=coef(summary(lm(biomass ~ metric, data = mcd_KNZ[ mcd_KNZ$bio_type==bt & mcd_KNZ$metric_type==mt, ])))[2,4]
    dats <- c(bt,mt,p.val)
    p.df <- rbind(p.df,dats)
  }
}
colnames(p.df) <- c("bio_type","metric_type","p.val")
p.df <- p.df[-1,]

mcd_KNZ_pval <- merge(mcd_KNZ,p.df)

## plot all diversity metrics and biomass
ggplot(data = mcd_KNZ_pval,
       aes(x=metric,y=biomass)) +
  facet_wrap(metric_type~bio_type, scales = "free", nrow = 7 , ncol = 5)+
  geom_smooth(data=subset(mcd_KNZ_pval,  p.val < 0.05), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, aes(color=as.numeric(p.val))) +
  geom_point() +
  theme_classic() +
  xlab("Diversity metric") +
  ylab("Biomass") +
  guides(color=guide_legend(title="p-value"))
ggsave(filename="Figures/Biomass_v_Diversity_KNZ.jpg",height=15,width=12)

###########################################################################################
###########################################################################################
##
## Combined plot
##
###########################################################################################
###########################################################################################

# add a site variable
mcd_SGS_pval$site <- rep("SGS",nrow(mcd_SGS_pval))
mcd_KNZ_pval$site <- rep("KNZ",nrow(mcd_KNZ_pval))

bothsites <- as.data.frame(rbind(mcd_SGS_pval,mcd_KNZ_pval))
# replace bogr/bohi or ange with 'dominant
bothsites$bio_type <- gsub("BOGR.BOHI","Dominant",bothsites$bio_type)
bothsites$bio_type <- gsub("ANGE","Dominant",bothsites$bio_type)
ggplot(data = bothsites,
       aes(x=metric,y=biomass)) +
  facet_wrap(metric_type~bio_type, scales = "free", nrow = 7 , ncol = 5)+
  geom_point(aes(shape=site)) +
  geom_smooth(data=subset(bothsites,  p.val < 0.05), method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, aes(color=as.numeric(p.val),lty=site)) +
  scale_shape_manual(values=c(1,16)) +
  theme_classic() +
  xlab("Diversity metric") +
  ylab("Biomass") +
  guides(color=guide_legend(title="p-value"))
ggsave(filename="Figures/Biomass_v_Diversity_bothsites.jpg",height=15,width=12)


###########################################################################################
###########################################################################################
##
## Illinois
##
###########################################################################################
###########################################################################################
# ## phylogenetic tree
# 
# community_ILL <- read.nexus("/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/phylogenetic-diversity-pilot/Illinois/Cipres_Data/RAxML_bestTree.result")
# community_ILL$edge.length <- round(community_ILL$edge.length, digits = 3)
# ggtree(community_ILL) +
#   geom_tiplab() +
#   geom_text(aes(label=node), hjust=-.3) +
#   geom_text(aes(label=branch.length), hjust=1, vjust=0.1) +
#   theme(plot.margin=unit(c(1,4,1,1),"cm")) # top right bottom left
# ggsave(width=25,height=15, filename="Illinois/Distances_toscale.jpg")
# 
# ###########################################################################################
# ## process species comp data
# 
# ILL_df <- ILL[,1:6]
# sort(unique(ILL_df$Species))#62 unique species 
# clean_ILL <- data.frame()
# for(i in 1:3){
#   for(j in 1:10){
#     PlotID <- ILL_df[(ILL_df$Plot == j),]
#     PlotID <- PlotID[(PlotID$Block == i),]
#     PlotID <- PlotID[order(PlotID[,'Species'],-PlotID[,'X.cover']),] ## sort so largest percent cover is on top
#     PlotID <- PlotID[!duplicated(PlotID$Species),] ## remove any species duplicated with lower percent cover
#     print(PlotID)
#     clean_ILL <- rbind(clean_ILL,PlotID)
#   }
# }
# 
# clean_ILL_wide <- clean_ILL[,-c(1,6)]
# clean_ILL_wide <- reshape(clean_ILL_wide, direction = "wide", timevar = "Species", idvar = c("Block","Plot"))
# clean_ILL_wide[is.na(clean_ILL_wide)] <- 0
# names(clean_ILL_wide) <- gsub("X.cover.", "", names(clean_ILL_wide))
# clean_ILL_wide <- clean_ILL_wide[,-c(1:2)]
# rownames(clean_ILL_wide) <- seq(1,20,1)
# 
# ## calculate Faith's PD (see picante package)
# tree_ILL <- prune.sample(clean_ILL_wide,community_ILL)
# phylo.div_ILL <- ses.pd(clean_ILL_wide,community_ILL, include.root = F)
# 
# ## Shannons + Simpsons's
# div_ILL_Shan <- community_diversity(clean_ILL, replicate.var = c("Plot","Block"), abundance.var = "X.cover", metric = "Shannon")
# div_ILL_Simp <- community_diversity(clean_ILL, replicate.var = c("Plot","Block"), abundance.var = "X.cover", metric = "InverseSimpson")
# 
# ## Berger Parker
# bp.index_ILL <- diversity(as.matrix(clean_ILL_wide), type = "berger-parker")
# bp.index_ILL <- bp.index_ILL[order(as.numeric(row.names(bp.index_ILL))),]
# 
# ## evenness  - would like to use a better metric..
# evenness_ILL_EQ <- community_structure(clean_ILL, replicate.var = c("Plot","Block"), abundance.var = "X.cover", metric = "EQ")
# evenness_ILL_Evar <- community_structure(clean_ILL, replicate.var = c("Plot","Block"), abundance.var = "X.cover", metric = "Evar")
# 
# ## concatenate
# total.div_ILL <- cbind(div_ILL_Shan,div_ILL_Simp,phylo.div_ILL,bp.index_ILL,evenness_ILL_EQ,evenness_ILL_Evar)
# write.csv(total.div_ILL, "DIVERSITY_ILL.csv")
# 
# ###########################################################################################
# ## import biomass data
# biomassdata = read.csv("biomass_compiled.csv",header=TRUE)
# biomassdata_ILL = biomassdata
# phylodata_ILL = read.csv("DIVERSITY_ILL.csv",header=TRUE)
# 
# nrow(biomassdata_ILL)
# nrow(phylodata_ILL)
# 
# completedata_ILL <-  cbind(biomassdata_ILL, phylodata_ILL)
# completedata_ILL <-  subset(completedata_ILL, select=-c(OLDLITTER,NEWLITTER,MOSS,TOTAL))
# completedata_ILL$Graminoids <-  completedata_ILL$INDGRASS + completedata_ILL$OTHRGRASS
# completedata_ILL$Total <- completedata_ILL$TOTALNOLITTER
# completedata_ILL$pd.obs.per.spp <- completedata_ILL$pd.obs / completedata_ILL$ntaxa ## calculate phylogenetic distance per species
# df.part_ILL <- subset(completedata_ILL, select=c(Block,Plot,INDGRASS,FORBS,Graminoids,Total,pd.obs.per.spp,richness,berger.parker.D,EQ,Evar,Shannon,InverseSimpson))
# ## weird outlier ..
# df.part_ILL <- df.part_ILL[-9,]
# mcd_ILL <- melt(df.part_ILL, id.vars = c("Block","Plot","pd.obs.per.spp","richness","berger.parker.D","EQ","Evar","Shannon","InverseSimpson"))
# mcd_ILL <- melt(mcd_ILL, id.vars = c("Block","Plot","variable","value"))
# colnames(mcd_ILL) <- c("Block","plot","bio_type","biomass","metric_type","metric")
# 
# ## plot all diversity metrics and biomass
# ggplot(data = mcd_ILL,
#        aes(x=metric,y=biomass)) +
#   facet_wrap(metric_type~bio_type, scales = "free", nrow = 7 , ncol = 4)+
#   geom_smooth(method=lm, se=F) +
#   geom_point() +
#   theme_classic() +
#   xlab("Diversity metric") +
#   ylab("Biomass")
# ggsave(filename="Figures/Biomass_v_Diversity_ILL.jpg",height=17,width=14)
# 
# ###########################################################################################
# ###########################################################################################
# ##
# ## Preliminary stats
# ##
# ###########################################################################################
# ###########################################################################################
# 
# summary(lm(Total~richness,completedata_KNZ))
# summary(lm(Total~pd.obs.per.spp,completedata_KNZ))
# summary(lm(Total~EQ,completedata_KNZ))
# summary(lm(Total~Evar,completedata_KNZ))
# summary(lm(Total~berger.parker.D,completedata_KNZ))
# summary(lm(Total~Shannon,completedata_KNZ))
# summary(lm(Total~InverseSimpson,completedata_KNZ))
# summary(lm(Total~richness + pd.obs.per.spp + EQ + berger.parker.D, completedata_KNZ))
# 

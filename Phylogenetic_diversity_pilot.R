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
  "/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity"
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

community_SGS <- read.nexus("/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/SGS/Cipres_Data/RAxML_bestTree.result")
community_SGS$edge.length <- round(community_SGS$edge.length, digits = 3)
#community.drop <- (subtrees(community)[[35]]) # don't want the outgroup
ggtree(community_SGS) +
  geom_tiplab() +
  geom_text(aes(label=node), hjust=-.3) +
  geom_text(aes(label=branch.length), hjust=1, vjust=0.1) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) # top right bottom left
ggsave(width=15,height=8, filename="SGS/Distances_toscale.jpg")

###########################################################################################
## process species comp data

SGS_df <- SGS[,1:5] # just keep the first five columns
SGS_df <- SGS_df[!(SGS_df$species == "Unk_6-_smallaster"),]
SGS_df <- SGS_df[!(SGS_df$species == "bare_ground"),]
SGS_df <- SGS_df[!(SGS_df$species == "cryptogramic"),]
SGS_df$species <- as.character(SGS_df$species)
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
clean_SGS_wide <- clean_SGS[,-1]
clean_SGS_wide <- reshape(clean_SGS_wide, direction = "wide", timevar = "species", idvar = c("Block","plot"))
clean_SGS_wide[is.na(clean_SGS_wide)] <- 0
names(clean_SGS_wide) <- gsub("X.cover.", "", names(clean_SGS_wide))
clean_SGS_wide <- clean_SGS_wide[,-c(1:2)]
rownames(clean_SGS_wide) <- seq(1,30,1)

## calculate Faith's PD (see picante package)
tree_SGS <- prune.sample(clean_SGS_wide,community_SGS)
phylo.div_SGS <- ses.pd(clean_SGS_wide,community_SGS, include.root = F)

## calculate other measures of diversity

## Shannons
clean_SGS$shan <- ( clean_SGS$X.cover / 100 ) * log( clean_SGS$X.cover / 100 )
shannon_SGS <- aggregate(shan ~ plot + Block, data=clean_SGS, FUN=sum)
shannon_SGS$shan <- -1 * shannon_SGS$shan

## Simpsons
clean_SGS$`n-1` <- clean_SGS$X.cover - 1
clean_SGS$`nxn-1` <- clean_SGS$`n-1` * clean_SGS$X.cover
all.simp_SGS <- aggregate(X.cover ~ Block + plot, data=clean_SGS, FUN=sum)
ind.simp_SGS <- aggregate(`nxn-1` ~ Block + plot, data=clean_SGS, FUN=sum)
simpson_SGS <- 1 - (ind.simp_SGS$`nxn-1` / (all.simp_SGS$X.cover * (all.simp_SGS$X.cover - 1)))

## Berger Parker
bp.index_SGS <- diversity(as.matrix(clean_SGS_wide), type = "berger-parker")
bp.index_SGS <- bp.index_SGS[order(as.numeric(row.names(bp.index_SGS))),]

## evenness  - would like to use a better metric..
evenness_SGS <- diversity(as.matrix(clean_SGS_wide, type="ev"))
evenness_SGS <- evenness_SGS[order(as.numeric(row.names(evenness_SGS))),]

## concatenate
total.div_SGS <- cbind(shannon_SGS,simpson_SGS,phylo.div_SGS,bp.index_SGS,evenness_SGS)
write.csv(total.div_SGS, "DIVERSITY_sgs.csv")

###########################################################################################
## import biomass data
biomassdata = read.csv("Biomass_data_all.csv",header=TRUE)
biomassdata_SGS = biomassdata[(biomassdata$Site == "SGS"),]
phylodata_SGS = read.csv("DIVERSITY_sgs.csv",header=TRUE)

nrow(biomassdata_SGS)
nrow(phylodata_SGS)

completedata_SGS <-  cbind(biomassdata_SGS, phylodata_SGS)
completedata_SGS <-  subset(completedata_SGS, select=-c(ANGE,Woody))
completedata_SGS$Graminoids <-  completedata_SGS$Grass
completedata_SGS$Total <-  completedata_SGS$total.cactus ## exclude cactus mass..
completedata_SGS$pd.obs.per.spp <- completedata_SGS$pd.obs / completedata_SGS$ntaxa ## calculate phylogenetic distance per species
df.part_SGS <- subset(completedata_SGS, select=c(Block,plot,BOGR.BOHI,Forbs,Graminoids,Total,pd.obs.per.spp,ntaxa,berger.parker.D,evenness,shan,simpson_SGS))
mcd_SGS <- melt(df.part_SGS, id.vars = c("Block","plot","pd.obs.per.spp","ntaxa","berger.parker.D","evenness","shan","simpson_SGS"))
mcd_SGS <- melt(mcd_SGS, id.vars = c("Block","plot","variable","value"))
colnames(mcd_SGS) <- c("Block","plot","bio_type","biomass","metric_type","metric")

## plot all diversity metrics and biomass
ggplot(data = mcd_SGS,
       aes(x=metric,y=biomass)) +
  facet_wrap(metric_type~bio_type, scales = "free", nrow = 6 , ncol = 4)+
  geom_smooth(method=lm, se=F) +
  geom_point() +
  theme_classic() +
  xlab("Diversity metric") +
  ylab("Biomass")
ggsave(filename="Figures/Biomass_v_Diversity_SGS.jpg",height=12,width=10)

###########################################################################################
###########################################################################################
##
## Konza
##
###########################################################################################
###########################################################################################
## phylogenetic tree

community_KNZ <- read.nexus("/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/Konza/Cipres_Data/RAxML_bestTree.result")
community_KNZ$edge.length <- round(community_KNZ$edge.length, digits = 3)
ggtree(community_KNZ) +
  geom_tiplab() +
  geom_text(aes(label=node), hjust=-.3) +
  geom_text(aes(label=branch.length), hjust=1, vjust=0.1) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) # top right bottom left
ggsave(width=25,height=15, filename="Konza/Distances_toscale.jpg")

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
phylo.div_KNZ <- ses.pd(clean_KNZ_wide,community_KNZ, include.root = F)

## calculate other measures of diversity

## Shannons
clean_KNZ$shan <- ( clean_KNZ$X.cover / 100 ) * log( clean_KNZ$X.cover / 100 )
shannon_KNZ <- aggregate(shan ~ plot + Block, data=clean_KNZ, FUN=sum)
shannon_KNZ$shan <- -1 * shannon_KNZ$shan

## Simpsons
clean_KNZ$`n-1` <- clean_KNZ$X.cover - 1
clean_KNZ$`nxn-1` <- clean_KNZ$`n-1` * clean_KNZ$X.cover
all.simp_KNZ <- aggregate(X.cover ~ Block + plot, data=clean_KNZ, FUN=sum)
ind.simp_KNZ <- aggregate(`nxn-1` ~ Block + plot, data=clean_KNZ, FUN=sum)
simpson_KNZ <- 1 - (ind.simp_KNZ$`nxn-1` / (all.simp_KNZ$X.cover * (all.simp_KNZ$X.cover - 1)))

## Berger Parker
bp.index_KNZ <- diversity(as.matrix(clean_KNZ_wide), type = "berger-parker")
bp.index_KNZ <- bp.index_KNZ[order(as.numeric(row.names(bp.index_KNZ))),]

## evenness  - would like to use a better metric..
evenness_KNZ <- diversity(as.matrix(clean_KNZ_wide, type="ev"))
evenness_KNZ <- evenness_KNZ[order(as.numeric(row.names(evenness_KNZ))),]

## concatenate
total.div_KNZ <- cbind(shannon_KNZ,simpson_KNZ,phylo.div_KNZ,bp.index_KNZ,evenness_KNZ)
write.csv(total.div_KNZ, "DIVERSITY_KNZ.csv")

###########################################################################################
## import biomass data
biomassdata = read.csv("Biomass_data_all.csv",header=TRUE)
biomassdata_KNZ = biomassdata[(biomassdata$Site == "KNZ"),]
phylodata_KNZ = read.csv("DIVERSITY_KNZ.csv",header=TRUE)

nrow(biomassdata_KNZ)
nrow(phylodata_KNZ)

completedata_KNZ <-  cbind(biomassdata_KNZ, phylodata_KNZ)
completedata_KNZ <-  subset(completedata_KNZ, select=-c(BOGR.BOHI,Cactus,total.cactus))
completedata_KNZ$Graminoids <-  completedata_KNZ$Grass + completedata_KNZ$ANGE
completedata_KNZ$Total <- completedata_KNZ$total
completedata_KNZ$pd.obs.per.spp <- completedata_KNZ$pd.obs / completedata_KNZ$ntaxa ## calculate phylogenetic distance per species
df.part_KNZ <- subset(completedata_KNZ, select=c(Block,plot,ANGE,Forbs,Graminoids,Woody,Total,pd.obs.per.spp,ntaxa,berger.parker.D,evenness,shan,simpson_KNZ))
mcd_KNZ <- melt(df.part_KNZ, id.vars = c("Block","plot","pd.obs.per.spp","ntaxa","berger.parker.D","evenness","shan","simpson_KNZ"))
mcd_KNZ <- melt(mcd_KNZ, id.vars = c("Block","plot","variable","value"))
colnames(mcd_KNZ) <- c("Block","plot","bio_type","biomass","metric_type","metric")

## plot all diversity metrics and biomass
ggplot(data = mcd_KNZ,
       aes(x=metric,y=biomass)) +
  facet_wrap(metric_type~bio_type, scales = "free", nrow = 6 , ncol = 5)+
  geom_smooth(method=lm, se=F) +
  geom_point() +
  theme_classic() +
  xlab("Diversity metric") +
  ylab("Biomass")
ggsave(filename="Figures/Biomass_v_Diversity_KNZ.jpg",height=15,width=14)

###########################################################################################
###########################################################################################
##
## Illinois
##
###########################################################################################
###########################################################################################
## phylogenetic tree

community_ILL <- read.nexus("/Users/avahoffman/Dropbox/Research/Genetic_Diversity_Pilot/Phylogenetic_diversity/Illinois/Cipres_Data/RAxML_bestTree.result")
community_ILL$edge.length <- round(community_ILL$edge.length, digits = 3)
ggtree(community_ILL) +
  geom_tiplab() +
  geom_text(aes(label=node), hjust=-.3) +
  geom_text(aes(label=branch.length), hjust=1, vjust=0.1) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) # top right bottom left
ggsave(width=25,height=15, filename="Illinois/Distances_toscale.jpg")

###########################################################################################
## process species comp data

ILL_df <- ILL[,1:6]
sort(unique(ILL_df$Species))#62 unique species 
clean_ILL <- data.frame()
for(i in 1:3){
  for(j in 1:10){
    PlotID <- ILL_df[(ILL_df$Plot == j),]
    PlotID <- PlotID[(PlotID$Block == i),]
    PlotID <- PlotID[order(PlotID[,'Species'],-PlotID[,'X.cover']),] ## sort so largest percent cover is on top
    PlotID <- PlotID[!duplicated(PlotID$Species),] ## remove any species duplicated with lower percent cover
    print(PlotID)
    clean_ILL <- rbind(clean_ILL,PlotID)
  }
}

clean_ILL_wide <- clean_ILL[,-c(1,6)]
clean_ILL_wide <- reshape(clean_ILL_wide, direction = "wide", timevar = "Species", idvar = c("Block","Plot"))
clean_ILL_wide[is.na(clean_ILL_wide)] <- 0
names(clean_ILL_wide) <- gsub("X.cover.", "", names(clean_ILL_wide))
clean_ILL_wide <- clean_ILL_wide[,-c(1:2)]
rownames(clean_ILL_wide) <- seq(1,20,1)

## calculate Faith's PD (see picante package)
tree_ILL <- prune.sample(clean_ILL_wide,community_ILL)
phylo.div_ILL <- ses.pd(clean_ILL_wide,community_ILL, include.root = F)

## calculate other measures of diversity

## Shannons
clean_ILL$shan <- ( clean_ILL$X.cover / 100 ) * log( clean_ILL$X.cover / 100 )
shannon_ILL <- aggregate(shan ~ Plot + Block, data=clean_ILL, FUN=sum)
shannon_ILL$shan <- -1 * shannon_ILL$shan

## Simpsons
clean_ILL$`n-1` <- clean_ILL$X.cover - 1
clean_ILL$`nxn-1` <- clean_ILL$`n-1` * clean_ILL$X.cover
all.simp_ILL <- aggregate(X.cover ~ Block + Plot, data=clean_ILL, FUN=sum)
ind.simp_ILL <- aggregate(`nxn-1` ~ Block + Plot, data=clean_ILL, FUN=sum)
simpson_ILL <- 1 - (ind.simp_ILL$`nxn-1` / (all.simp_ILL$X.cover * (all.simp_ILL$X.cover - 1)))

## Berger Parker
bp.index_ILL <- diversity(as.matrix(clean_ILL_wide), type = "berger-parker")
bp.index_ILL <- bp.index_ILL[order(as.numeric(row.names(bp.index_ILL))),]

## evenness  - would like to use a better metric..
evenness_ILL <- diversity(as.matrix(clean_ILL_wide, type="ev"))
evenness_ILL <- evenness_ILL[order(as.numeric(row.names(evenness_ILL))),]

## concatenate
total.div_ILL <- cbind(shannon_ILL,simpson_ILL,phylo.div_ILL,bp.index_ILL,evenness_ILL)
write.csv(total.div_ILL, "DIVERSITY_ILL.csv")

###########################################################################################
## import biomass data
biomassdata = read.csv("biomass_compiled.csv",header=TRUE)
biomassdata_ILL = biomassdata
phylodata_ILL = read.csv("DIVERSITY_ILL.csv",header=TRUE)

nrow(biomassdata_ILL)
nrow(phylodata_ILL)

completedata_ILL <-  cbind(biomassdata_ILL, phylodata_ILL)
completedata_ILL <-  subset(completedata_ILL, select=-c(OLDLITTER,NEWLITTER,MOSS,TOTAL))
completedata_ILL$Graminoids <-  completedata_ILL$INDGRASS + completedata_ILL$OTHRGRASS
completedata_ILL$Total <- completedata_ILL$TOTALNOLITTER
completedata_ILL$pd.obs.per.spp <- completedata_ILL$pd.obs / completedata_ILL$ntaxa ## calculate phylogenetic distance per species
df.part_ILL <- subset(completedata_ILL, select=c(Block,Plot,INDGRASS,FORBS,Graminoids,Total,pd.obs.per.spp,ntaxa,berger.parker.D,evenness,shan,simpson_ILL))
## weird outlier ..
df.part_ILL <- df.part_ILL[-9,]
mcd_ILL <- melt(df.part_ILL, id.vars = c("Block","Plot","pd.obs.per.spp","ntaxa","berger.parker.D","evenness","shan","simpson_ILL"))
mcd_ILL <- melt(mcd_ILL, id.vars = c("Block","Plot","variable","value"))
colnames(mcd_ILL) <- c("Block","plot","bio_type","biomass","metric_type","metric")

## plot all diversity metrics and biomass
ggplot(data = mcd_ILL,
       aes(x=metric,y=biomass)) +
  facet_wrap(metric_type~bio_type, scales = "free", nrow = 6 , ncol = 4)+
  geom_smooth(method=lm, se=F) +
  geom_point() +
  theme_classic() +
  xlab("Diversity metric") +
  ylab("Biomass")
ggsave(filename="Figures/Biomass_v_Diversity_ILL.jpg",height=15,width=14)


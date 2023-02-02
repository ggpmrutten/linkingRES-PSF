## Collects root traits from the database based on species list
##  also adds phylogenetic tree and a figure
##
##  Project: Review Plant & Soil
##
##  by  Gemma Rutten (gemma.rutten@unibe.ch)
##      Last Edited January 23
##
##
### Used data from Petermann et al 2008 and Bennett et al 

## clean working space 
cat("\014") 
rm(list=ls())

## load packages
library(vegan)
library(ggplot2)
library(grid)
library(gridExtra)
library(ape)
library(geiger)
library(tidyverse)
library(phytools)

## functions
mansca<-function(x){x/sqrt(sum(x^2)/(length(x)-1))}

## Read sRoot modified root traits for which we have phylogeny
load("Data/sRoot.traits.Rdata") #cleaned according to Bergmann et al from sRoot traits
summary(spectrait)# 568 species 
str(spectrait)

## Read pruned tree for those species 
tree.prp <- read.tree("Data/sRoot_species.tre")
length(tree.prp$tip.label) # 568 species

# See which species don't match original species list and phylogeny tips
setdiff(tree.prp$tip.label, spectrait$full_species) # should be character(0)
setdiff(spectrait$full_species, tree.prp$tip.label) # should be character(0)

## Read species lists and PSF data
species.RES.pet <-read.csv2("Data/linking.RES.PSF.SpeciesList.pet.csv")# Petermann 24
species.RES.pet$full_species <- gsub(" ", "_", species.RES.pet$species.full)
species.RES.ben<-read.csv2("Data/linking.RES.PSF.SpeciesList.ben.csv")# Bennett 44
species.RES.ben$full_species <- gsub(" ", "_", species.RES.ben$species.full)

#spectrait[grep("Melilotus", spectrait$full_species), c("full_species")]
# bind species lists together
species.RES<-rbind(species.RES.pet, species.RES.ben)#68

### now start links to traits, selecting complete observations and links to phylogeny
# merge with trait data
traits.RES<-merge(species.RES, spectrait, by="full_species", all.x=T) 
summary(traits.RES)# 68 species, 33 NA's

# remove species without trait data drop NAs 
traits.RES<-traits.RES %>% 
  drop_na("Root_diameter_corrected")%>%
  drop_na("rootN_corrected")%>%
  drop_na("RTD_corrected")%>%
  drop_na("SRL_corrected")
str(traits.RES)## for 35 species these traits are in the database

## subset phylo tree
# matrix with species names as rows
species.names <- as.matrix(cbind(
  full_species=traits.RES$full_species                        , 
  species.full=traits.RES$species.full))

dimnames(species.names)[[1]] <- species.names[,1]

str(species.names)##35

## check if names in tree and data are same  
comb<-geiger::name.check(tree.prp,species.names)

## remove tips from tree
tree.RES<-drop.tip(tree.prp, comb$tree_not_data)
comb$data_not_tree

geiger::name.check(tree.RES,species.names)
length(tree.RES$tip.label)#35
str(species.names)#35

## save pruned tree
write.tree(tree.RES,file="Data/treeRES.tre")

## phylogenetically corrected PCA
## make sure names are the same
setdiff(tree.RES$tip.label, traits.RES$full_species) # should be character(0)
setdiff(traits.RES$full_species, tree.RES$tip.label) # should be character(0)
rownames(traits.RES) <- traits.RES$full_species

# PCA with phy tree
phylpca <- phytools::phyl.pca(tree.RES, traits.RES[,4:7], mode="corr", method="lambda")
summary(phylpca)
print(phylpca)

# save PCA axes loadings
traits.RES$PCA1<-phylpca$S[,1];traits.RES$PCA1rev<--phylpca$S[,1];traits.RES$PCA2<-phylpca$S[,2]##

## clean levels of growthform
traits.RES$growthForm<- gsub("shrub/tree", "tree", traits.RES$growthForm)
traits.RES[traits.RES$nitrogenFixationNodDB=="YES",]$growthForm<- gsub("herb", "legume", traits.RES[traits.RES$nitrogenFixationNodDB=="YES",]$growthForm)
traits.RES$micro<-as.factor(paste(traits.RES$growthForm,traits.RES$mycorrhizalAssociationTypeFungalRoot, sep="_"))

table(as.factor(traits.RES$growthForm))
table(as.factor(traits.RES$micro))# less traits for AMF trees than EMF and 4,4,5 for g,h,l
## save data file
save (traits.RES, file = "Data/traits.RES.benpet.Rdata")## all together 

## correlations
names(traits.RES)
summary(lm(PCA1 ~ Root_diameter_corrected, data=traits.RES))#R2 adj 0.8116, p-value: 1.013e-13
summary(lm(PCA2 ~ Root_diameter_corrected, data=traits.RES))#R2 adj 0.2917 , p-value: 0.0004814

summary(lm(PCA1 ~ SRL_corrected, data=traits.RES))#R2 adj 0.9369, p-value: < 2.2e-16
summary(lm(PCA2 ~ SRL_corrected, data=traits.RES))#R2 adj  0.1, p-value: 0.03599

summary(lm(PCA1 ~ rootN_corrected, data=traits.RES))#R2 adj -0.03021, p-value: NS
summary(lm(PCA2 ~ rootN_corrected, data=traits.RES))#R2 adj 0.4588, p-value: 4.724e-06

summary(lm(PCA1 ~ RTD_corrected, data=traits.RES))#R2 adj -0.01814 , p-value: NS
summary(lm(PCA2 ~ RTD_corrected, data=traits.RES))#R2 adj 0.463 , p-value: 4.144e-06

######  Figure: All Species in RES  #####
gemma<-c(16,17)
gems<-c("olivedrab3", "lightgoldenrod4")
summary(phylpca)

## first calculate vectors 
# fit $ loadings is fit1 $ rotation is phy.pca$L
# fit $ scores is fit1 $ x is phy.pca$S
datapc <- data.frame(varnames=rownames(phylpca$L), phylpca$L)#rotation
mult <- min(
  (max(phylpca$S[,2]) - min(phylpca$S[,2])/(max(phylpca$L[,2])-min(phylpca$L[,2]))),
  (max(phylpca$S[,1]) - min(phylpca$S[,1])/(max(phylpca$L[,1])-min(phylpca$L[,1])))
)
datapc1 <- transform(datapc,
                     v1 = .7 * mult * (get("PC1")),
                     v2 = .7 * mult * (get("PC2")))
datapc1$lab<-c("D","SRL","RTD","N")

# save figure 
tiff('Plots/Together.in.traitspaceCorrected.tiff',
     width = 20, height = 20, units = "cm", res = 400 , pointsize= 15,bg="transparent")
par(mar = c(4, 4, 1, 1),xpd = T)

#par(mar=c(3,3,3,3))
#plot(-phylpca$S[order(row.names(phylpca$S)),1], phylpca$S[order(row.names(phylpca$S)),2], pch=19, col="gray",
#     xlab=paste0("PC1"," (", round(0.4913097 ,2),")"), ylab=paste0("PC2"," (", round(0.2680714 ,2),")"),
#     main=NA,  xlim=c(-50,50),ylim=c(-50,50) ,cex=0.4)

#points(phylpca$S[which(species.means.4$mycorrhizalAssociationTypeFungalRoot=="ErM"),1],
#       -phylpca$S[which(species.means.4$mycorrhizalAssociationTypeFungalRoot=="ErM"),2],
#       pch=20, col="grey80",cex=1)

plot(traits.RES$PCA1rev,traits.RES$PCA2, pch=gemma[droplevels(as.factor(traits.RES$mycorrhizalAssociationTypeFungalRoot))], 
     col=gems[droplevels(as.factor(traits.RES$woodiness))],
     xlim=c(-50,50),ylim=c(-50,50),
     xlab=paste0("PC1"," (", round(0.4913097 ,2),")"),
     ylab=paste0("PC2"," (", round(0.2680714 ,2),")"))
#abline(h=0, col="darkgray",xpd = F); abline(v=0,col="darkgray",xpd = F)
polygon(c(53.65,-53.65,-53.65,53.65),c(0,0,53.65,53.65),density=30,  angle=45, col= rgb(255,0,0,max=255,alpha = (100 - 80) * 255 / 100))#pathogen
polygon(c(-53.65,-53.65,0,0),c(53.65,-53.65,-53.65,53.65),density=30, angle=-45,col= rgb(0,0,255,max=255,alpha = (100 - 80) * 255 / 100))#mutualists
text(traits.RES$PCA1rev,traits.RES$PCA2, font=3, col="black",
     labels = traits.RES$species.full,pos=1, offset = 0.15, cex=.5)

arrows(0,0,-datapc1$v1[1],datapc1$v2[1],col="black",length = 0.05, angle = 35)
arrows(0,0,-datapc1$v1[2],datapc1$v2[2],col="black",length = 0.05, angle = 35)
arrows(0,0,-datapc1$v1[3],datapc1$v2[3],col="black",length = 0.05, angle = 35)
arrows(0,0,-datapc1$v1[4],datapc1$v2[4],col="black",length = 0.05, angle = 35)

text(-datapc1$v1,datapc1$v2,pos=c(2,4,1,3),offset =.5, font=1, labels = datapc1$lab,  cex=.8, col="black") 

text(-40,33,'outsourcer-fast', col=rgb(1,0,1), cex=1)
text(40,33,'DIY-fast', col=rgb(1,0,0),cex=1)
text(-40,-33,'outsourcer-slow',col=rgb(0,0,1), cex=1)
text(40,-33, 'DIY-slow', col=c("darkgray"),cex=1)

legend("topright",inset = c( -0, 0),  
       paste0(levels(droplevels(as.factor(traits.RES$woodiness)))," (",table(traits.RES$woodiness), ")"),
       cex=1, pch=16,col=gems, bty="n" )
legend("bottomright",inset = c( -0, 0),
       paste0(levels(droplevels(as.factor(traits.RES$mycorrhizalAssociationTypeFungalRoot))),
              " (",table(droplevels(as.factor(traits.RES$mycorrhizalAssociationTypeFungalRoot))), ")"),
       cex=1, pch=gemma, bty="n")
dev.off()

# save figure 
tiff('Plots/concept.traitspace.tiff',
     width = 20, height = 20, units = "cm", res = 400 , pointsize= 15,bg="transparent")
par(mar = c(3, 3, 1, 1),xpd = T)

plot(traits.RES$PCA1rev,traits.RES$PCA2, col="white",
     xlim=c(-50,50),ylim=c(-50,50),
     xlab="", ylab="")
abline(h=0, col="darkgray",xpd = F); abline(v=0,col="darkgray",xpd = F)
polygon(c(53.65,-53.65,-53.65,53.65),c(0,0,53.65,53.65),density=30,  angle=45, col= rgb(255,0,0,max=255,alpha = (100 - 80) * 255 / 100))#pathogen
polygon(c(-53.65,-53.65,0,0),c(53.65,-53.65,-53.65,53.65),density=30, angle=-45,col= rgb(0,0,255,max=255,alpha = (100 - 80) * 255 / 100))#mutualists
#text(-traits.RES$PCA1,traits.RES$PCA2, font=3,
#     labels = traits.RES$species.full,pos=1, offset = 0.1, cex=.5)

arrows(0,0,-datapc1$v1[1],datapc1$v2[1],col="black",length = 0.05, angle = 35, lwd=2)
arrows(0,0,-datapc1$v1[2],datapc1$v2[2],col="black",length = 0.05, angle = 35, lwd=2)
arrows(0,0,-datapc1$v1[3],datapc1$v2[3],col="black",length = 0.05, angle = 35, lwd=2)
arrows(0,0,-datapc1$v1[4],datapc1$v2[4],col="black",length = 0.05, angle = 35, lwd=2)

text(-datapc1$v1,datapc1$v2,pos=c(2,4,1,3),offset =.5, font=1, labels = datapc1$lab,  cex=.8, col="black") 

text(-40,33,'outsourcer-fast', col=rgb(1,0,1), cex=1)
text(40,33,'DIY-fast', col=rgb(1,0,0),cex=1)
text(-40,-33,'outsourcer-slow',col=rgb(0, 0, 1), cex=1)
text(40,-33, 'DIY-slow', col=c("darkgray"),cex=1)

dev.off()
# fast fig to check independent of reversed order of the PCA1

p.pca<-phyl.pca(tree.RES,traits.RES[,4:7],method="lambda", mode="corr")
summary(p.pca)
biplot(p.pca)

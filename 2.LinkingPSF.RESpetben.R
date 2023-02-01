## Combines root traits with PSF measures
## and calculates pairwise differences of PSFs in Root Economics spectrum 
##
##  Review Plant&Soil
##
##  by  Gemma Rutten (gemma.rutten@unibe.ch)
##      Last Edited January 23
##
##
###

## clean working space 
cat("\014") 
rm(list=ls())

## set path
#path <- "C:/analysis/linkingRES&PSF"
#setwd(path)

## load packages
library(vegan)
library(ggplot2)
library(grid)
library(gridExtra)
library(ape)
#library(geiger)
library(tidyverse)

## load functions 
len<-function(x){length(na.omit(x))} 
mansca<-function(x){x/sqrt(sum(x^2)/(length(x)-1))}

## tweaks 

## load data

# Table with plant species names and performance on home and away soil

#load("Data/home.vs.Away.PSF.cra.Rdata") #PSFs crawford
#summary(PSFs.A)# 776 comparisons
#names(PSFs.A)# Crawford has reciprocal for A and B= PSFpaiwise
#"BinB.mean","BinB.se","BinB.N" , "BinA.mean", "BinA.se", "BinA.N", "rrIs",
#"Variance.rrIs." ,  "HAb","pwAB","var"
#sel.PSFs.A<-PSFs.A %>% 
#  select("Species.pair","Species.A", "Family.A", "Species.B", "Family.B","AinA.mean", "AinA.se",
#         "AinA.N", "AinB.mean", "AinB.se","AinB.N","HAa",)
#sel.PSFs.A$Study<-c("crawford")## selected common columns

# Table with plant species names and performance on home and away soil
load("Data/home.vs.Away.PSF.ben.Rdata") #PSFs bennett
summary(PSFs.B)# 79 comparisons
names(PSFs.B)# Bennett has sterile treatment for cons and hetero separately =  PSFlife/control
# "st.AinA.mean", "st.AinA.se", "st.AinA.N",    "st.AinB.mean", "st.AinB.se",   "st.AinB.N","Hbio","Abio"
sel.PSFs.B<-PSFs.B %>% 
  select(-c("Family.A"  ,   "mycorrhiza.A", "Family.B"  ,   "mycorrhiza.B"))
sel.PSFs.B$Study<-c("bennett")## selected common columns

# Table with plant species names and performance on home and away soil
# save(PSFs.C, file = "Data/home.vs.Away.PSF.pet.Rdata") 
load("Data/home.vs.Away.PSF.pet.Rdata") #PSFs petermann
summary(PSFs.C)# 240 comparisons with 5 different treatments
names(PSFs.C)# Peterman has a control and sterile treatment

lifePSFs.C<-PSFs.C %>% 
  filter(treat == "control")%>%
  select("Species.pair","Species.A", "Species.B", "AinA.mean", "AinA.se",
         "AinA.N", "AinB.mean", "AinB.se","AinB.N","HAa")

sterPSFs.C<-PSFs.C %>% 
  filter(treat == "sterile" )%>%
  select("Species.pair","AinA.mean", "AinA.se",
         "AinA.N", "AinB.mean", "AinB.se","AinB.N")%>%
  `colnames<-`(c("Species.pair", "st.AinA.mean", "st.AinA.se",
                 "st.AinA.N", "st.AinB.mean", "st.AinB.se","st.AinB.N"))

#save(Pairwise.Petermann, file = "Data/pw.PSF.pet.Rdata") #pairwise
load("Data/pw.PSF.pet.Rdata") #PSFs petermann
summary(Pairwise.Petermann)# 240 comparisons with 5 different treatments
names(Pairwise.Petermann)# Peterman has a control and sterile treatment

lifePSFs.C<-Pairwise.Petermann %>% 
  filter(treat == "control")%>%
  select("Species.pair","Species.A", "Species.B", "AinA.mean", "AinA.se",
         "AinA.N", "AinB.mean", "AinB.se","AinB.N","BinB.mean","BinB.se","BinB.N"
         , "BinA.mean", "BinA.se", "BinA.N","HAa", "HAb", "pairwise")

sterPSFs.C<-Pairwise.Petermann %>% 
  filter(treat == "sterile" )%>%
  select("Species.pair","AinA.mean", "AinA.se",
         "AinA.N", "AinB.mean", "AinB.se","AinB.N")%>%
  `colnames<-`(c("Species.pair", "st.AinA.mean", "st.AinA.se",
                 "st.AinA.N", "st.AinB.mean", "st.AinB.se","st.AinB.N"))

sel.PSFs.C<-merge(lifePSFs.C,sterPSFs.C, by="Species.pair")

sel.PSFs.C$Hbio = log(sel.PSFs.C$AinA.mean / sel.PSFs.C$st.AinA.mean)
sel.PSFs.C$Abio = log(sel.PSFs.C$AinB.mean / sel.PSFs.C$st.AinB.mean)
sel.PSFs.C$Study<-c("petermann")## selected common columns

# combine data sets bennett and petermann
# NOTE: we leave A out because there is no sterile treatment
#names(sel.PSFs.A)# crawford
names(sel.PSFs.B);length(sel.PSFs.B$Species.pair)#79
names(sel.PSFs.C[,-c(10:15,17,18)]);length(sel.PSFs.C[,-c(10:15,17,18)]$Species.pair)#48

PSFs<-rbind(sel.PSFs.B,sel.PSFs.C[,-c(10:15,17,18)])
str(PSFs)#127

# Tables with Species and selected root traits
## save (traits.RES, file = paste0(root, "/traits.RES.benpet.Rdata")) in sRoot_extr RES script
#load( "C:/analysis/GRooT-Data-master/traits.RES.benpet.Rdata")#
#load( "C:/analysis/GRooT-Data-master/DataFiles/RES.Species.together.Rdata")
load( "Data/traits.RES.benpet.Rdata")#
summary(traits.RES)# 35 species with root traits in the database
names(traits.RES)

## pairwise pearson correlations between 4 traits
# source("C:/analysis/functions/correlationplots.R")
# sdat<-traits.RES%>%
#  select(Root_diameter_corrected,RTD_corrected,rootN_corrected,SRL_corrected)
# r_plot<-round(cor(sdat,use="pairwise.complete.obs",method= "pearson"),2)
# write.table(r_plot, file="corr.coef.csv", sep=",")

# trial<-round(cor.test.p(sdat), 3)
# write.table(trial, file="corr.p.csv", sep=",")
#  D ~ SRL seems ok r=-0.53 p= <0.001
#  N ~ RTD not at all correlated

## bind together bind PSF and Traits
try2<-PSFs %>% 
  mutate(species.full=Species.A) %>% 
  left_join(traits.RES, by="species.full") %>% 
  mutate(species.full=Species.B) %>% 
  left_join(traits.RES, by="species.full") %>%   
  drop_na("SRL_corrected.x")%>%
  drop_na("SRL_corrected.y")%>%
  drop_na("rootN_corrected.x")%>%
  drop_na("rootN_corrected.y")

names(try2)
summary(try2)## 61 pairwise comparisons
head(try2)

### Calculate PSF direction based root traits using these formula
## PSF    <- ifelse(away>home,(home/away-1),(home+abs(away)/abs(away)))
## PSFlog <- ifelse(away>home,log(home/away),log(home+abs(away)/abs(away)))
## PSF    <- away – home / ABS(home)
try3 <- try2 %>%
  mutate(coll.diff = (PCA1rev.y-PCA1rev.x),
         cons.diff = (PCA2.y-PCA2.x),
         coll.mid = ((try2$PCA1rev.x + try2$PCA1rev.y)/2),
         cons.mid = ((PCA2.x + PCA2.y)/2),
         #pairwise.RES = coll_grad + cons_grad,
         #
         cons_fact = as.factor(ifelse(PCA2.x < 0 , "slow","fast")),
         coll_fact = as.factor(ifelse(PCA1rev.x < 0 ,"outsource","DIY")),
         #
         cons_trait = as.factor(ifelse(cons.mid < 0 , "slow","fast")),
         coll_trait = as.factor(ifelse(coll.mid < 0 ,"outsource","DIY")),
         strategytrait = as.factor(paste(coll_trait,cons_trait,sep="-")),
         #
         myco.stat = as.factor(paste(mycorrhizalAssociationTypeFungalRoot.x,mycorrhizalAssociationTypeFungalRoot.y,sep="_")),
         GFm.stat = as.factor(paste(micro.x,micro.y,sep="_")),
         form.stat = as.factor(paste(growthForm.x,growthForm.y,sep="_")), 
         strategy = as.factor(paste(coll_fact,cons_fact,sep="_")))
str(try3)# 61 obs.
#sd(try3$cons.diff)/sqrt(length(try3$cons.diff))

#add scaled variables
#try3$sc.coll.mid <- scale (try3$coll.mid)
#try3$sc.coll.diff <- scale (try3$coll.diff)
#try3$sc.cons.mid <- scale (try3$cons.mid)
#try3$sc.cons.diff <- scale (try3$cons.diff)

## save clean combined file
save (try3, file ="Data/PSF.RES.benpet.Rdata") # 2.linking PSF with RES

## save a species list for supplementary
str(try3)
tableS2<-try3 %>%
  select(Species.pair,coll.diff,cons.diff,coll.mid,cons.mid,HAa)
write.table(tableS2, "tableS2.csv")

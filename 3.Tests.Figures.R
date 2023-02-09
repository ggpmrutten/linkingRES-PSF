## Stats and figures Linking RES and PSF
##  Review Plant&Soil
##
##  by  Gemma Rutten (gemma.rutten@unibe.ch)
##      Last Edited February 23
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
library(ggeffects) 
library(lme4)
library(RColorBrewer)

## load functions 
# source("C:/analysis/functions/GreenFigures.R")
len<-function(x){length(na.omit(x))} 
mansca<-function(x){x/sqrt(sum(x^2)/(length(x)-1))}

## load tweaks 
myArrow <- function(x0, y0, x1, y1, cut = .95, ...){
  x.new <- (1 - cut) * x0 + cut * x1
  y.new <- (1 - cut) * y0 + cut * y1
  # segments(x0, y0, x1, y1, ...) #if you want the arrow to continue
  arrows(x0, y0, x.new, y.new, ...)
}

## set colors
wGem<-c("#FFFFE5")

yGem<-c("#FE9929")
oGem<-c("#CC4C02")
rGem<-c("#662506")# col2rgb(rGem)

gemma<-c(16,17)
gems<-c("#35978F","#807DBA")#colors for symbols woody non-woody

## load data
# save (try3, file = paste0(path, "/Data/PSF.RES.benpet.Rdata")) # 2.linking PSF with RES
load("Data/PSF.RES.benpet.Rdata") #load combined data
names(try3)
table(try3$Study)
##

## test linear correlations between the variables
#source("C:/analysis/functions/correlationplots.R")
#sdat<-try3%>%
#  select(SRL_corrected.x,SRL_corrected.y,rootN_corrected.x,rootN_corrected.y,PCA1rev.x, PCA2.x,PCA1rev.y, PCA2.y,coll.mid,coll.diff,cons.mid,cons.diff)
#r_plot<-round(cor(sdat,use="pairwise.complete.obs",method= "pearson"),2)
#write.table(r_plot, file="corr.coef.derived.vars.csv", sep=",")

## best fit for correlations linear or quadratic?
y<-try3$coll.mid
x<-try3$coll.diff; x2<-x^2
lin   <-lm(y~ x)
quad  <-lm(y~ x+ x2)### calculate quad model I(x^2)

anova(lin, quad)# no linear corr but quad fit for conservation only.(pairwise also for)

## stat tests
## HA effects explained by species pairs characteristics

# full model Feedback Home.away
m1<-lm(HAa~ 
       + coll.mid
       + cons.mid
       + coll.diff
       + cons.diff
       + coll.mid : cons.mid
       + coll.mid : coll.diff
       + coll.mid : cons.diff
       + cons.mid : coll.diff
       + cons.mid : cons.diff
       + coll.diff: cons.diff
       + coll.mid : cons.mid : coll.diff
       + coll.mid : cons.mid : cons.diff 
       + cons.mid : coll.diff: cons.diff
       + coll.mid : coll.diff: cons.diff
       + coll.mid : cons.mid : coll.diff: cons.diff# -6.136e-06  1.836e-06  -3.343  0.00168 **
       , data=try3)

summary(m1)
anova(m1)
# sign 4 way interaction
# Residual standard error: 0.4486 on 45 degrees of freedom
# Multiple R-squared:  0.5721,	Adjusted R-squared:  0.4295 
# F-statistic: 4.011 on 15 and 45 DF,  p-value: 0.0001479

# full model Feedback Home sterile
m2<-lm(Hbio~ 
       + PCA1rev.x
       + PCA2.x              # 11.4588 0.001294 **
       + PCA1rev.x :PCA2.x   #  5.9355 0.017986 * 
       , data=try3)
anova(m2)
summary(m2)
# Residual standard error: 0.7028 on 57 degrees of freedom
# Multiple R-squared:  0.2348,	Adjusted R-squared:  0.1946 
# F-statistic: 5.831 on 3 and 57 DF,  p-value: 0.001514

## figures

# PSFs in Trait space (Figure 4)
tiff('Plots/PSFs.in.traitspace.petben.tiff',
     width = 20, height = 20, units = "cm", res = 400 , pointsize= 15,bg="white")
par(mar = c(4, 4, 1, 1),xpd = T)

#try3$Homecols<-ifelse(try3$Hbio<0, "red", "green")

plot(try3$PCA1rev.x,try3$PCA2.x,  
     col=gems[droplevels(as.factor(try3$woodiness.x))],
     pch=16,
     xlim=c(-50,50),ylim=c(-50,50),
     xlab=c("Collaboration gradient"),
     ylab=c("Conservation gradient")
)
#segments(try3$PCA1rev.x,try3$PCA2.x,-try3$PCA1.y,try3$PCA2.y, col=c("lightgray"))
polygon(c(53.65,-53.65,-53.65,53.65),c(0,0,53.65,53.65), density=30, angle=45, col= rgb(102,37,6, max=255, alpha = 51))#pathogen
polygon(c(-53.65,-53.65,0,0),c(53.65,-53.65,-53.65,53.65),density=30, angle=-45, col= rgb(254, 153, 41, max = 255, alpha = 51))#mutualists

abline(h=0, col="lightgray",xpd = F)
abline(v=0, col="lightgray",xpd = F)

myArrow(try3$PCA1rev.x,try3$PCA2.x,try3$PCA1rev.y,try3$PCA2.y, length = 0.1, angle = 25,
       col=c("darkgray"))

text(try3$PCA1rev.x,try3$PCA2.x, font=3, 
     labels = try3$Species.A, pos=1, offset = .15, cex=0.5, col="darkgray") 

segments(-expl.cons$PCA1.x,expl.cons$PCA2.x,-expl.cons$PCA1.y,expl.cons$PCA2.y, col=c("green"))
#segments(-expl.coll$PCA1.x,expl.coll$PCA2.x,-expl.coll$PCA1.y,expl.coll$PCA2.y, col=c("orange"))

text(-40,33,'outsourcer-fast', col=oGem, cex=1)
text(40,33,'DIY-fast', col=rGem,cex=1)
text(-40,-33,'outsourcer-slow',col=yGem, cex=1)
text(40,-33, 'DIY-slow', col=c("darkgray"),cex=1)

points(try3$PCA1rev.y,try3$PCA2.y,  
       col=gems[droplevels(as.factor(try3$woodiness.y))],pch=1)

#legend("topright",inset = c( 0, 0),  
#       paste0(levels((as.factor(try3$strategy)))," (",table(as.factor(try3$strategy)), ")"),
#       col=stratcols, cex=1, pch=16, bty="n" )
legend("topright",inset = c( -0, 0),  
       paste0(levels(droplevels(as.factor(try$woodiness)))," (",table(try$woodiness), ")"),
       cex=1, pch=16,col=gems, bty="n" )

legend("bottomright",inset = c(0, 0),c("Home","Away"),pch=c(16,1),bty="n")

points(try3$PCA1rev.x,try3$PCA2.x,  
       col=gems[droplevels(as.factor(try3$woodiness.x))],pch=1)

points(try3$PCA1rev.x,try3$PCA2.x,  
       col=gems[droplevels(as.factor(try3$woodiness.x))],pch=16)

dev.off()

# relations between variables (supplementary)
stratcols<-c(rGem, "gray", oGem,yGem)

tiff('Plots/relationsbetweenvariables.tiff',
     width = 30, height = 15, units = "cm", res = 400 , pointsize= 12,bg="white")
if (interactive()) {
  split.screen(c(1,2))       # split display into 2 screens
  
  screen(1) # prepare screen 1 for output##Adjust for area
  par(oma=c(0,0,0,0),mar=c(4,4,1,1))
  plot(cons.mid~cons.diff , data=try3, pch=16, las=1,bty="l", tck=.005, 
       ylab="conservation gradient (location)", xlab="conservation gradient (distance)",cex.lab=1,
       col=stratcols[strategy],col.lab="black",col.axis="black")
  
  y<-try3$cons.mid 
  x<-try3$cons.diff
  x2<-x^2
  quad  <-lm(y~ x+ x2)### calculate quad model I(x^2)
  
  mq<-summary(quad)
  R2.quad<-round(mq$adj.r.squared,digit=2)                       ### substract R2quad
  P.quad<-round(pf(mq$fstatistic[1], mq$fstatistic[2], mq$fstatistic[3], lower.tail = F),digit=3) ### substract Pvalues
  ax<-summary(x); ay<-summary(y)
  text(ax[6]*0.75,ay[6]*0.95,substitute(paste("R"^2, " = ", R2.quad,"***", sep=" "),list(R2.quad=R2.quad)), col="black", cex=1.2)
  xv<-seq(ax[1],ax[6], 0.1)## for quadratic
  yv<-predict(quad,list(x=xv,x2=xv^2))## for quadratic
  lines(xv,yv,lwd=1, col="black")        ## smooth for quadratic
  
  
  screen(2) # prepare screen 2 for output
  par(oma=c(0,0,0,0),mar=c(4,4,1,1))
  plot(coll.mid ~ coll.diff, data=try3, pch=16, las=1,bty="l", tck=.005, 
       ylab="collaboration gradient (location)", xlab="collaboration gradient (distance)",cex.lab=1,
       col=stratcols[strategy],col.lab="black",col.axis="black")
  
  y<-try3$coll.mid
  x<-try3$coll.diff
  lin  <-lm(y~ x)#
  mq<-summary(lin)
  R2.quad<-abs(round(mq$adj.r.squared,digit=2))                       ### substract R2quad
  P.quad<-round(pf(mq$fstatistic[1], mq$fstatistic[2], mq$fstatistic[3], lower.tail = F),digit=3) ### substract Pvalues
  ax<-summary(x); ay<-summary(y)
  text(ax[6]*0.75,ay[6]*0.95,substitute(paste("R"^2, " = ", R2.quad," n.s.", sep=" "),list(R2.quad=R2.quad)), col="black", cex=1.2)
  
  close.screen(all = TRUE) 
}
dev.off()

# PSF home sterile over gradients (Figure 5)
table(try3$Species.A)
# several measures per home so first average
meanHbio<-try3 %>%
  group_by(Species.A)%>%
  summarize(PCA1=mean(PCA1rev.x),PCA2=mean(PCA2.x),PSF=mean(Hbio))
write.table(meanHbio, "tableS1.csv")## suppl S1          

# categorize Feedback strength            
range(meanHbio$PSF)
range(try3$Hbio)

meanHbio$PSF_fact<-cut(meanHbio$PSF, c(-3,-2,-1,0,1), labels= c(-2,-1,0,1))

# # PSF home sterile over gradients (Figure 5) the actual figure
psfcols<-c("darkred", "red", "lightgray","darkblue")

tiff('Plots/IsolatedHomeEffects.tiff',
     width = 15, height = 15, units = "cm", res = 400 , pointsize= 12,bg="white")
par(oma=c(0,0,0,0),mar=c(4,4,1,1))

plot(meanHbio$PCA1, meanHbio$PCA2, col=psfcols[meanHbio$PSF_fact],pch=16,
     xlim=c(-50,50),ylim=c(-50,50),
     xlab=c("Collaboration gradient"),
     ylab=c("Conservation gradient"))
abline(v=0, col="lightgray");abline(h=0, col="lightgray")

text(-try3$PCA1.x,try3$PCA2.x, font=3,
     labels = try3$Species.A ,pos=3, offset = 0.15, cex=.5, col="darkgray")

legend("topright",inset = c( -0, 0), title=expression("PSF"[live/control]),
       paste0(levels(droplevels(as.factor(meanHbio$PSF_fact)))),
       cex=1, pch=16,col=psfcols, bty="n" )

text(-30,30,'outsourcer-fast', col=oGem, cex=1)
text(30,30,'DIY-fast', col=rGem,cex=1)
text(-30,-30,'outsourcer-slow',col=yGem, cex=1)
text(30,-30, 'DIY-slow', col=c("darkgray"),cex=1)

dev.off()

# show the three way interactions (Figure 6)

## now separate the continuous variables into 3 categories
try3$cons.diff_fact<-as.factor(cut(try3$cons.diff,c(-55,-5,5,55), labels=c("negative", "same strategy", "positive")))
table(try3$cons.diff_fact)

try3$coll.diff_fact<-as.factor(cut(try3$coll.diff,c(-55,-5,5,55), labels=c("negative", "same strategy", "positive")))
table(try3$coll.diff_fact)

## the figure
tiff('Plots/4wayInteraction_x_colldistRange.tiff',
     width = 45, height = 15, units = "cm", res = 400 , pointsize= 22,bg="white")
if (interactive()) {
  split.screen(c(1,3)) # split display into 2 screens
  
  screen(1) # prepare screen 1 for output##Adjust for area
  par(oma=c(0,0,0,0),mar=c(4,5,1,1))
  
  # no distance cons gradient
  plot(HAa~coll.diff,data=try3[try3$cons.diff_fact=="same strategy",],
       col=stratcols[try3[try3$cons.diff_fact=="same strategy",]$strategytrait], 
       main="a) Same on conservation axis", pch=16, ylim=c(-2.5,1),las=1,bty="l", tck=.005, 
       xlab="", ylab= expression("PSF"[home/away]),cex.lab=1,cex.main = 1, 
       font.main=1 ,col.lab="black",col.axis="black")
  abline(a=0,b=0, col="lightgray")
  #anova(lm(HAa~coll.diff.srl,data=try3[try3$cons.diff_fact=="same strategy",]))
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="same strategy",])
  xs <- range(try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="same strategy",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rGem, lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="same strategy",])
  xs <- range(try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="same strategy",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = "gray", lwd = 2)
  
  #mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="same strategy",])
  #xs <- range(try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="same strategy",]$coll.diff)
  #ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  #lines(xs, ys, col = "#66CD00", lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="same strategy",])
  xs <- range(try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="same strategy",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = yGem, lwd = 2)
  
  screen(2) # prepare screen 1 for output##Adjust for area 
  par(oma=c(0,0,0,0),mar=c(4,3,1,1))
  
  # negative distance cons gradient
  plot(HAa~coll.diff, col=stratcols[try3[try3$cons.diff_fact=="negative",]$strategytrait], data=try3[try3$cons.diff_fact=="negative",],
       main="b) Negative distance", pch=16 ,ylim=c(-2.5,1),las=1,bty="l", tck=.005, 
       xlab="Collaboration distance", ylab="",cex.lab=1,cex.main = 1, font.main=1 ,col.lab="black",col.axis="black")
  abline(a=0,b=0, col="lightgray")
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rGem, lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = "gray", lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = oGem, lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = yGem, lwd = 2)
  
  screen(3) # prepare screen 1 for output##Adjust for area
  par(oma=c(0,0,0,0),mar=c(4,3,1,1))
  # positive distance cons gradient
  plot(HAa~coll.diff, col=stratcols[try3[try3$cons.diff_fact=="positive",]$strategytrait], data=try3[try3$cons.diff_fact=="positive",], 
       main="c) Positive distance", pch=16 ,ylim=c(-2.5,1),las=1,bty="l", tck=.005, 
       xlab="", ylab="",cex.lab=1,col.lab="black",cex.main = 1, font.main=1 ,col.axis="black")
  abline(a=0,b=0, col="lightgray")
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="positive",])
  xs <- range(try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="positive",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rGem, lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="positive",])
  xs <- range(try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="positive",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = "gray", lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="positive",])
  xs <- range(try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="positive",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = oGem, lwd = 2)
  
  #mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="positive",])
  #xs <- range(try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="positive",]$coll.diff)
  #ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  #lines(xs, ys, col = "darkgreen", lwd = 2)
  
  close.screen(all = TRUE) 
}
dev.off()
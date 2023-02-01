## Stats and figures Linking RES and PSF
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
library(ggeffects) 
library(lme4)

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
gems<-c("olivedrab3", "lightgoldenrod4")

# PSFs in Trait space (Figure 4)
tiff('Plots/PSFs.in.traitspace.petben.tiff',
     width = 20, height = 20, units = "cm", res = 400 , pointsize= 15,bg="white")
par(mar = c(4, 4, 1, 2),xpd = T)

#try3$Homecols<-ifelse(try3$Hbio<0, "red", "green")

plot(try3$PCA1rev.x,try3$PCA2.x,  
     col=gems[droplevels(as.factor(try3$woodiness.x))],
     pch=16,
     xlim=c(-50,50),ylim=c(-50,50),
     xlab=c("Collaboration gradient"),
     ylab=c("Conservation gradient")
)
#segments(try3$PCA1rev.x,try3$PCA2.x,-try3$PCA1.y,try3$PCA2.y, col=c("lightgray"))

abline(h=0, col="lightgray",xpd = F)
abline(v=0, col="lightgray",xpd = F)

myArrow(try3$PCA1rev.x,try3$PCA2.x,try3$PCA1rev.y,try3$PCA2.y, length = 0.1, angle = 25,
       col=c("darkgray"))

text(try3$PCA1rev.x,try3$PCA2.x, font=3, 
     labels = try3$Species.A, pos=1, offset = .15, cex=0.5, col="darkgray") 

#segments(-expl.cons$PCA1.x,expl.cons$PCA2.x,-expl.cons$PCA1.y,expl.cons$PCA2.y, col=c("green"))
#segments(-expl.coll$PCA1.x,expl.coll$PCA2.x,-expl.coll$PCA1.y,expl.coll$PCA2.y, col=c("orange"))

text(-40,33,'outsourcer-fast', col=rgb(1,0,1), cex=1)
text(40,33,'DIY-fast', col=rgb(1,0,0),cex=1)
text(-40,-33,'outsourcer-slow',col=rgb(0, 0, 1), cex=1)
text(40,-33, 'DIY-slow', col=c("darkgray"),cex=1)

points(try3$PCA1rev.y,try3$PCA2.y,  
       col=gems[droplevels(as.factor(try3$woodiness.y))],pch=1)

#legend("topright",inset = c( 0, 0),  
#       paste0(levels((as.factor(try3$strategy)))," (",table(as.factor(try3$strategy)), ")"),
#       col=stratcols, cex=1, pch=16, bty="n" )
legend("bottomright",inset = c(0, 0),c("Home","Away"),pch=c(16,1),bty="n")

points(try3$PCA1rev.x,try3$PCA2.x,  
       col=gems[droplevels(as.factor(try3$woodiness.x))],pch=1)

points(try3$PCA1rev.x,try3$PCA2.x,  
       col=gems[droplevels(as.factor(try3$woodiness.x))],pch=16)

dev.off()

# relations between variables (supllementary)
stratcols<-c(rgb(1,0, 0), "gray", rgb(1, 0, 1),rgb(0, 0, 1))

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

# the actual figure
psfcols<-c("red", "orange", "lightgray","darkblue")

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

text(-40,33,'outsourcer-fast', col=rgb(1,0,1), cex=1)
text(40,23,'DIY-fast', col=rgb(1,0,0),cex=1)
text(-40,-33,'outsourcer-slow',col=rgb(0, 0, 1), cex=1)
text(40,-33, 'DIY-slow', col=c("darkgray"),cex=1)

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
  plot(HAa~coll.diff,data=try3[try3$cons.diff_fact=="same strategy",],col=stratcols[try3[try3$cons.diff_fact=="same strategy",]$strategytrait], 
       main="a) Same on conservation axis", pch=16, ylim=c(-2.5,1),las=1,bty="l", tck=.005, 
       xlab="", ylab= expression("PSF"[home/away]),cex.lab=1,cex.main = 1, 
       font.main=1 ,col.lab="black",col.axis="black")
  abline(a=0,b=0, col="lightgray")
  #anova(lm(HAa~coll.diff.srl,data=try3[try3$cons.diff_fact=="same strategy",]))
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="same strategy",])
  xs <- range(try3[try3$strategytrait=="DIY-fast" & try3$cons.diff_fact=="same strategy",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rgb(1,0,0), lwd = 2)
  
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
  lines(xs, ys, col = rgb(0,0,1), lwd = 2)
  
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
  lines(xs, ys, col = rgb(1,0,0), lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = "gray", lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rgb(1,0,1), lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="negative",])
  xs <- range(try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="negative",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rgb(0,0,1), lwd = 2)
  
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
  lines(xs, ys, col = rgb(1,0,0), lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="positive",])
  xs <- range(try3[try3$strategytrait=="DIY-slow" & try3$cons.diff_fact=="positive",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = "gray", lwd = 2)
  
  mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="positive",])
  xs <- range(try3[try3$strategytrait=="outsource-fast" & try3$cons.diff_fact=="positive",]$coll.diff)
  ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  lines(xs, ys, col = rgb(1,0,1), lwd = 2)
  
  #mO <- lm(HAa~coll.diff,data=try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="positive",])
  #xs <- range(try3[try3$strategytrait=="outsource-slow" & try3$cons.diff_fact=="positive",]$coll.diff)
  #ys <- predict(mO, newdata = data.frame(coll.diff = xs))
  #lines(xs, ys, col = "darkgreen", lwd = 2)
  
  close.screen(all = TRUE) 
}
dev.off()

## figure with fitted values for factors
##predict values with unscaled data results are the same 
m.fac<-lm(HAa~ 
          + coll_fact
          + cons_fact
          + coll.diff
          + cons.diff_fact
          + coll_fact : cons_fact
          + coll_fact : coll.diff
          + coll_fact : cons.diff_fact
          + cons_fact : coll.diff
          + cons_fact : cons.diff_fact
          + coll.diff: cons.diff_fact
          + coll_fact : cons_fact :coll.diff
          + coll_fact : cons_fact :cons.diff_fact #1 90.619 9.604 0.001942 **
          + cons_fact : coll.diff :cons.diff_fact
          + coll_fact : coll.diff :cons.diff_fact
          + coll_fact : cons_fact: coll.diff :cons.diff_fact
          , data=try3)
names(try3)
anova(m.fac)

m.fac<-lm(HAa~ 
          + strategytrait
          + coll.diff
          + cons.diff_fact
          + strategytrait : coll.diff
          + strategytrait : cons.diff_fact
          + coll.diff: cons.diff_fact
          + strategytrait: coll.diff :cons.diff_fact
          , data=try3)
names(try3)
anova(m.fac)


## predict with ggpredict
library(ggeffects)
tiff('Plots/4wayinteractionwithggeffects.jpeg', width = 15, height = 9, units = "cm", res = 400 , pointsize= 12,bg="white")
ggpredict(m.fac, c("coll.diff", "strategytrait","cons.diff_fact")) |> plot(add.data = T, limit.range = T)+
  labs( x = c("Collaboration distance"), y = c("PSF (home/away)"), title="PSF by root strategy", colour="conservation axis (midpoint)")
dev.off()

  labs(
    colour = get_legend_title(dat),
    x = get_x_title(dat),
    y = get_y_title(dat),
    title = get_title(dat)
  )


ggplot(gg.prd, aes(coll.mid, predicted))
  geom_point(position = position_dodge(.1)) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    position = position_dodge(.1)
  ) +
  theme_bw()+
  xlab("delta conservation gradient")+ylab("Plant soil Feedback (Home Away)")+
  #coord_cartesian(xlim=c(0,2.5),ylim=c(0,4))+
  theme(legend.position="bottom")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = rel(1.1)),
        panel.border = element_blank(),axis.line = element_line()) 
  scale_x_discrete(breaks = 1:3, labels = get_x_labels(dat))

ggplot(d, aes(x=cons.diff,y= DV))+
  geom_point(data=try3,aes(x=cons.diff,y=scale(HAa)))+
  #geom_point(size=1, shape=1, color ="black")+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  xlab("delta conservation gradient")+ylab("Plant soil Feedback (Home Away)")+
  #coord_cartesian(xlim=c(0,2.5),ylim=c(0,4))+
  theme(legend.position="bottom")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = rel(1.1)),
        panel.border = element_blank(),axis.line = element_line()) 

d = with(try3, 
         data.frame(coll.mid = rep(seq(min(coll.mid), max(coll.mid), length.out=20),6),
                    cons.mid = rep(seq(min(cons.mid), max(cons.mid), length=5),each=24),
                    cons.diff= rep(seq(min(cons.diff),max(cons.diff), length=3),20)))
d$fitHAa = predict(m2, newdata=d)

# add factors
d$cons_fact <- as.factor(ifelse(d$cons.mid < 0 , "slow", "fast"))
d$coll_fact <- as.factor(ifelse(d$coll.mid < 0 ,"mycorrhiza","root"))
d$cons.diff_fact <- as.factor(cut(d$cons.diff,3))
levels(d$cons.diff_fact)<-c("small","middle","large")

gems<-c("red", "gray")
gemma<-c(1,16,17)
par(oma=c(0,0,0,0),mar=c(4,4,1,1))
plot(fitHAa ~ coll.mid, data=d, pch=gemma[d$cons.diff_fact], las=1,bty="l", tck=.005, 
     xlab="collaboration gradient (mid point)", ylab="fitted values of feedback (HAa)",cex.lab=1,
     col=gems[d$cons_fact],col.lab="black",col.axis="black")

abline(a=0,b=0, col="lightgray")

abline(lm(fitHAa~coll.mid,data=d[d$cons_fact=="slow" & d$cons.diff_fact=="small",]), col="lightgray",lwd=2)
abline(lm(fitHAa~coll.mid,data=d[d$cons_fact=="slow" & d$cons.diff_fact=="middle",]), col="snow4",lwd=2)
abline(lm(fitHAa~coll.mid,data=d[d$cons_fact=="slow" & d$cons.diff_fact=="large",]), col="black",lwd=2)

abline(lm(fitHAa~coll.mid,data=d[d$cons_fact=="fast" & d$cons.diff_fact=="small",]), col="orange",lwd=2)
abline(lm(fitHAa~coll.mid,data=d[d$cons_fact=="fast" & d$cons.diff_fact=="middle",]),col="orangered",lwd=2)
abline(lm(fitHAa~coll.mid,data=d[d$cons_fact=="fast" & d$cons.diff_fact=="large",]), col="red4",lwd=2)

abline(coefs[1], coefs[2],
       col = "red",
       lwd = 1.5)

plot(cons.mid ~ cons.diff, data=d, col=gems[d$cons_fact],pch=gemma[d$cons.diff_fact], data=d)
abline(a=0,b=0, col="lightgray")
abline(lm(resi.PSF~resi.coll.mid,data=try3[try3$cons_fact=="fast",]), col="red")
abline(lm(resi.PSF~resi.coll.mid,data=try3[try3$cons_fact=="slow",]), col="gray")

gems1<-c("darkgreen", "green")

#plot data
figa<-ggplot(d, aes(x=cons.diff,y= DV))+
  geom_point(data=try3,aes(x=cons.diff,y=scale(HAa)))+
  #geom_point(size=1, shape=1, color ="black")+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  xlab("delta conservation gradient")+ylab("Plant soil Feedback (Home Away)")+
  #coord_cartesian(xlim=c(0,2.5),ylim=c(0,4))+
  theme(legend.position="bottom")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = rel(1.1)),
        panel.border = element_blank(),axis.line = element_line()) 
figa

## fancier test accounting for SpeciesA while no reps for species pairs
library(lme4)
#library(lmerTest)
## fancier test accounting for Species A, while no reps for species pairs

#library(lme4)
#source("C:/Analysis/functions/stepba function.R")

## HA effects explained by species pairs characteristics
m1<-lmer(scale(HAa)~ 
           + coll.mid
         + cons.mid
         + cons.diff
         + coll.diff
         + coll.mid:cons.mid
         + coll.mid:cons.diff
         + coll.mid:coll.diff
         + cons.mid:cons.diff
         + cons.mid:coll.diff
         + cons.diff:coll.diff
         + coll.mid:cons.mid:coll.diff
         + coll.mid:cons.mid:cons.diff #1 90.619 9.604 0.001942 **
         + cons.mid:cons.diff:coll.diff
         + (1|Species.A)
         , data=try3)
head(try3)
source("C:/Analysis/functions/stepba function.R")
stepba(m1)      # significant 3way interaction

m2<-lmer(HAa~ 
           + woodiness.x
         + coll.mid#0.01159 *
         #+ cons.mid
         + cons.diff
         + coll.diff
         #+ woodiness.x:coll.mid# 0.07946 .
         + woodiness.x:coll.diff
         #+ woodiness.x:cons.mid
         + woodiness.x:cons.diff#   0.0636 .
         #+ coll.mid:cons.mid
         + cons.diff:coll.diff
         + coll.mid:coll.diff
         #+ cons.mid:cons.diff
         + (1|Species.A)
         , data=try3)

anova(m1,m2)

qqnorm(residuals(m1)); qqline(residuals(m1)) 

# forests
anova(lm(HAa~coll.mid*cons.mid*cons.diff*coll.diff, data=try3[try3$woodiness.x=="woody",]))# 37 PSFs forests

# grasslands
anova(lm(HAa~coll.mid*cons.mid*cons.diff*coll.diff, data=try3[try3$woodiness.x=="non-woody",]))# 24 PSFs grasslands
anova(lm(HAa~coll.mid*cons.mid*coll.diff, data=try3[try3$woodiness.x=="non-woody",]))# 24 PSFs grasslands


#plot(lm(scale(HAa)~woodiness.x*coll.mid*cons.mid*cons.diff, data=try3))#
m1<-lm(HAa~coll.mid, data=try3)# correct PSFs for coll mid data 
try3$resid.m1<-residuals(m1)

m2<-lm(resid.m1~woodiness.x*cons.diff, data=try3)
anova(m2)

# create new data with all predictors like in model, the variable of interest gets range like in raw data
# all other continuous predictors are set to mean for categories include all levels
newdat <- expand.grid(cons.diff=seq(-2.525369,  2.525369,length.out=100),
                      coll.mid=mean(try3$coll.mid),
                      woodiness.x=c(levels(as.factor(try3$woodiness.x))))

#predict data                      
newdat$fit  <- predict(m1,newdat,allow.new.levels=T, na.action=na.pass, re.form=NA)# exclude random factors

gems1<-c("darkgreen", "green")
#plot data
figa<-ggplot(try3, aes(x=cons.diff,y= resid.m1))+
  geom_point(data=try3,aes(x=cons.diff,y=resid.m1),col=gems1[as.factor(try3$woodiness.x)],size=1)+
  #geom_point(size=1, shape=1, color ="black")+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  xlab("delta conservation gradient")+ylab("Plant soil Feedback (Home Away)")+
  #coord_cartesian(xlim=c(0,2.5),ylim=c(0,4))+
  theme(legend.position="bottom")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = rel(1.1)),
        panel.border = element_blank(),axis.line = element_line()) 
figa

#pars model

# create new data with all predictors like in model, the variable of interest gets range like in raw data
# all other continuous predictors are set to mean for categories include all levels
newdat <- expand.grid(cons.diff=seq(-2.525369,  2.525369,length.out=100),
                      coll.mid=mean(try3$coll.mid),
                      woodiness.x=c(levels(as.factor(try3$woodiness.x))))

#predict data                      
newdat$fit  <- predict(m1,newdat,allow.new.levels=T, na.action=na.pass, re.form=NA)# exclude random factors

gems1<-c("darkgreen", "green")
plot(try3$Hbio~try3$cons.diff, col=gems1[as.factor(try3$woodiness.x)], data=try3)
abline(a=0,b=0)
abline(lm(HAa~cons.diff, data=try3[try3$woodiness.x=="woody",]), col="green")
abline(lm(HAa~cons.diff, data=try3[try3$woodiness.x=="non-woody",]), col="darkgreen")


anova(lm(HAa~coll.mid*cons.mid*cons.diff*coll.diff, data=try3[try3$woodiness.x=="woody",]))# 37 PSFs forests
plot(lm(HAa~coll.mid*cons.mid*cons.diff, data=try3[try3$woodiness.x=="woody",]))

# Analysis of Variance Table
# Response: HAa
#                               Df Sum Sq Mean Sq F value    Pr(>F)    
#  coll.mid                     1 4.8878  4.8878 15.3781 0.0004949 ***
#  cons.mid                     1 0.0719  0.0719  0.2262 0.6379474    
#  cons.diff                    1 2.4423  2.4423  7.6840 0.0096297 ** 
#  coll.mid:cons.mid            1 0.1003  0.1003  0.3156 0.5785910    
#  coll.mid:cons.diff           1 1.0588  1.0588  3.3312 0.0782929 .  
#  cons.mid:cons.diff           1 0.0168  0.0168  0.0528 0.8199112    
#  coll.mid:cons.mid:cons.diff  1 1.8049  1.8049  5.6786 0.0239442 *  
#  Residuals                   29 9.2175  0.3178                      
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(lm(HAa~coll.mid*cons.mid*coll.diff, data=try3[try3$woodiness.x=="non-woody",]))# 24 PSFs grasslands
#plot(lm(HAa~coll.mid*cons.mid*coll.diff, data=try3[try3$woodiness.x=="non-woody",]))# 24 PSFs grasslands

# Analysis of Variance Table
# Response: HAa
#                              Df Sum Sq Mean Sq F value    Pr(>F)    
#  coll.mid                     1 0.42580 0.42580 12.9915 0.002377 **
#  cons.mid                     1 0.19310 0.19310  5.8918 0.027387 * 
#  coll.diff                    1 0.17675 0.17675  5.3929 0.033734 * 
#  coll.mid:cons.mid            1 0.00162 0.00162  0.0494 0.826987   
#  coll.mid:coll.diff           1 0.00270 0.00270  0.0823 0.777899   
#  cons.mid:coll.diff           1 0.00059 0.00059  0.0180 0.894899   
#  coll.mid:cons.mid:coll.diff  1 0.23400 0.23400  7.1395 0.016702 * 
#  Residuals                   16 0.52440 0.03277                    
#  ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## subset grassland species
grass<-try3[try3$woodiness.x=="non-woody",]
grass$grassLL<-as.factor(paste(grass$growthForm.x,grass$growthForm.y, sep="_"))
anova(lm(HAa~coll.mid*cons.mid*coll.diff, data=grass))# 24 PSFs grasslands

## take out variation explained by conservation gradient
g1<-lm(HAa~cons.mid, data=grass)# correct PSFs for cons mid data 
grass$resid.g1<-residuals(g1)

m2<-lm(resid.g1~coll.mid*coll.diff, data=grass)
anova(m2)
#plot(m2)

gems<-c(1:3)
plot(resid.g1~coll.diff, col=gems[as.factor(grass$growthForm.y)], data=grass)
plot(resid.g1~coll.mid, col=gems[as.factor(grass$growthForm.x)],data=grass)
head(grass)
abline(a=0,b=0)

m2<-lm(resid.m1~woodiness.x*cons.diff, data=try3)
anova(m2)

# create new data with all predictors like in model, the variable of interest gets range like in raw data
# all other continuous predictors are set to mean for categories include all levels
newdat <- expand.grid(cons.diff=seq(-2.525369,  2.525369,length.out=100),
                      coll.mid=mean(try3$coll.mid),
                      woodiness.x=c(levels(as.factor(try3$woodiness.x))))

#predict data                      
newdat$fit  <- predict(m1,newdat,allow.new.levels=T, na.action=na.pass, re.form=NA)# exclude random factors

library(ggplot2)
gems1<-c("darkgreen", "green")
#plot data
figa<-ggplot(try3, aes(x=cons.diff,y= resid.m1))+
  geom_point(data=try3,aes(x=cons.diff,y=resid.m1),col=gems1[as.factor(try3$woodiness.x)],size=1)+
  #geom_point(size=1, shape=1, color ="black")+
  geom_smooth(method = "lm",color="black")+
  theme_bw()+
  xlab("delta conservation gradient")+ylab("Plant soil Feedback (Home Away)")+
  #coord_cartesian(xlim=c(0,2.5),ylim=c(0,4))+
  theme(legend.position="bottom")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = rel(1.1)),
        panel.border = element_blank(),axis.line = element_line()) 
figa

## test correlations between the derived variables
sdat<-try3%>%
  select(coll.mid,coll.diff,cons.mid,cons.diff )
r_plot<-round(cor(sdat,use="pairwise.complete.obs",method= "pearson"),2)
write.table(r_plot, file="corr.coef.derived.vars.csv", sep=",")
round(cor.test.p(as.matrix(sdat)), 3)

m2<-aov(Hbio~strategy.x*micro.x, data=try3)#more balanced
summary(m2)
TukeyHSD(m2,"strategy.x")
plot(try3$Hbio~try3$strategy.x)
abline(a=0,b=0)

m3<-aov(Abio~strategy.y, data=try3[try3$strategy.x=="root.slow",] )#more balanced
summary(m3)
TukeyHSD(m3,"strategy.x")
gems<-c(1,2)
plot(try3$Abio~try3$strategy.x)
abline(a=0,b=0)

names(try3)
## very unbalanced to non woody AMF not nitFix so we select those
try4<-try3 %>% 
  filter(myco.statII=="AM_AM"& wood.stat== "non-woody"& Nfix.stat=="none")%>% 
  droplevels()
table(try4$form.statII)

## easy straight forward test
names(try4)

## Home away effects explained by coll or con diff or diff in mycotype
m1<-lm(HAa~ 
         #+ form.statII
         + cons.diff
       + coll.diff
       + cons.diff:coll.diff
       , data=try4)

anova(m1)# does not look good!

## fancier test accounting for SpeciesA while no reps for species pairs
#library(lmerTest)
## HA effects explained by species pairs characteristics
m1<-lmer(HAa~ 
           #+ form.statII
           + cons.diff
         + coll.diff
         #+ form.statII:cons.diff
         #+ form.statII:coll.diff
         + cons.diff:coll.diff
         #+ form.statII:cons.diff:coll.diff
         + (1|Species.A)
         + (1|Species.B)
         , data=try4)

m2<-lmer(HAa~ 
           + form.statII
         #+ cons.diff
         + coll.diff
         #+ form.statII:cons.diff
         + form.statII:coll.diff
         #+ cons.diff:coll.diff
         #+ form.statII:cons.diff:coll.diff
         + (1|Experiment)
         #+ (1|Study) no difference
         , data=try4)

qqnorm(residuals(m1)); qqline(residuals(m1)) 
anova(m1,m2)
#source("C:/Analysis/functions/stepba function.R")
#stepba(m1)      # significant 3way interaction

anova(m1)

# Type III Analysis of Variance Table with Satterthwaite's method
#                        Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# form.statII           2.40180 0.80060     3 124.68  4.3331 0.006091 **
# cons.diff             1.34400 1.34400     1 124.97  7.2741 0.007961 **
# coll.diff             1.40988 1.40988     1 125.10  7.6307 0.006604 **
# form.statII:cons.diff 1.59481 0.53160     3 118.19  2.8772 0.039025 * 
# form.statII:coll.diff 1.64693 0.54898     3 124.12  2.9712 0.034442 * 
# cons.diff:coll.diff   0.12647 0.12647     1 124.93  0.6845 0.409627   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(m2)

## the graphs
try3$micro.x

# tweaks
g4col=c("green", "darkgreen", "yellow", "orange")
g5col=c("yellow", "orange", "red","green","darkgreen")
g2col=c("green","darkgreen")
g2pch<-c(1,16)

# subset examples for both gradients
expl.coll<-try3 %>%
  filter(abs(cons.diff)<.3)

expl.cons<-try3 %>%
  filter(abs(coll.diff)<.3)
# traits instead of PCA and PSFs instead of species
tiff('Plots/PSFs.in.traitspace.petben.tiff',
     width = 20, height = 20, units = "cm", res = 400 , pointsize= 15,bg="white")
par(mar = c(4, 4, 1, 2),xpd = T)

plot(try3$SRL_corrected.x,try3$rootN_corrected.x, pch=c(16), 
     col=g5col[(as.factor(try3$micro.x))],
     xlim=c(-2,2),ylim=c(-2,2),
     xlab=c("Specific Root Length"),
     ylab=c("Nitrogen content")
)
abline(h=0, col="darkgray",xpd = F)
abline(v=0, col="darkgray",xpd = F)

segments(try3$SRL_corrected.x,try3$rootN_corrected.x,try3$SRL_corrected.y,try3$rootN_corrected.y, col=c("lightgray"))
segments(expl.coll$SRL_corrected.x,expl.coll$rootN_corrected.x,expl.coll$SRL_corrected.y,expl.coll$rootN_corrected.y, col=c("orange"))
segments(expl.cons$SRL_corrected.x,expl.cons$rootN_corrected.x,expl.cons$SRL_corrected.y,expl.cons$rootN_corrected.y, col=c("green"))

text(-1.5,1.5,'outsourcer-fast', col=c("darkgreen"), cex=.8)
text(1.5,1.5,'DIY-fast', col=c("orange"),cex=.8)
text(-1.5,-1.5,'outsourcer-slow',col=c("green"), cex=.8)
text(1.5,-1.5, 'DIY-slow', col=c("darkgray"),cex=.8)

points(try3$SRL_corrected.y,try3$rootN_corrected.y,  col=g5col[(as.factor(try3$micro.y))], pch=c(1))
legend("topright",inset = c( 0, 0),  
       paste0(levels(as.factor(try3$micro.x))," (",table(as.factor(try3$micro.x)), ")"),
       col=g5col, cex=1, pch=16, bty="n" )
legend("bottomright",inset = c(0, 0),c("Species A","Species B"),pch=c(16,1),bty="n")
dev.off()

## Does cons gradient explain PSF? YES!
spec.col=c("green", "darkgreen", "yellow", "orange")
plot(try3$HAa ~ try3$cons.diff , pch=16,col=spec.col[as.factor(try3$myco.stat)])
anova(lm(HAa~cons.diff*myco.stat, data=try3))
abline(lm(HAa ~cons.diff, data=try3))
abline(h=0, col="gray")


## Does coll gradient explain PSF? NO but EcM show positive and AM negative feedbacks
## even when they do not share mycorrhiza
plot(HAa ~ coll.diff, data=try3  , pch=16,col=spec.col[as.factor(try3$myco.stat)])
anova(lm(HAa~coll.diff*myco.stat, data=try3))
m3<-aov(HAa ~ coll.diff*myco.stat, data=try3)
boxplot(HAa ~ myco.stat, data=try3,col=spec.col)
TukeyHSD(m3, "myco.stat")
abline(h=0, col="gray")


## a figure: traits instead of PCA and PSFs instead of species
###a figure with PCA axes RES
expl.coll<-try3 %>%
  filter(abs(cons.diff)<1)

expl.cons<-try3 %>%
  filter(abs(coll.diff)<1)

gems<-c("green","darkgreen")
tiff('Plots/PSFs.in.traitspace.ben.tiff',
     width = 20, height = 20, units = "cm", res = 400 , pointsize= 15,bg="white")
par(mar = c(4, 4, 1, 2),xpd = T)

plot(try3$sc.SRL.x,try3$sc.N.x, pch=c(16), 
     col=gems[(as.factor(try3$mycosimple.x))],
     xlim=c(-3,3),ylim=c(-3,3),
     xlab=c("Specific Root Length"),
     ylab=c("Nitrogen content"))
abline(h=0, col="darkgray",xpd = F)
abline(v=0, col="darkgray",xpd = F)

segments(try3$sc.SRL.x,try3$sc.N.x,try3$sc.SRL.y,try3$sc.N.y, col=c("lightgray"))
segments(expl.coll$sc.SRL.x,expl.coll$sc.N.x,expl.coll$sc.SRL.y,expl.coll$sc.N.y, col=c("orange"))
segments(expl.cons$sc.SRL.x,expl.cons$sc.N.x,expl.cons$sc.SRL.y,expl.cons$sc.N.y, col=c("green"))

text(-2,2,'outsourcer-fast', col=c("darkgreen"), cex=.8)
text(2,2,'DIY-fast', col=c("orange"),cex=.8)
text(-2,-2,'outsourcer-slow',col=c("green"), cex=.8)
text(2,-2, 'DIY-slow', col=c("darkgray"),cex=.8)

points(try3$sc.SRL.y,try3$sc.N.y,  col=gems[(as.factor(try3$mycosimple.y))], pch=c(1))
legend("topright",inset = c( 0, 0),  
       paste0(levels((as.factor(try3$mycosimple.x)))," (",table(try3$mycosimple.x), ")"),
       col=gems, cex=1, pch=16, bty="n" )
legend("bottomright",inset = c( 0, 0),c("Species A","Species B"),pch=c(16,1),bty="n")
dev.off()

## Does cons gradient explain PSF? YES!expl.cons
spec.col=c("green", "darkgreen", "yellow", "orange")
plot(HAa ~ cons.diff ,data=expl.cons, pch=16,col=spec.col[as.factor(expl.cons$myco.statII)])
anova(lm(HAa~cons.diff*myco.statII, data=expl.cons))
abline(lm(HAa ~cons.diff, data=expl.cons))
abline(h=0, col="gray")

## Does coll gradient explain PSF? expl.coll
plot(HAa ~ coll.diff, data=expl.coll  , pch=16,col=spec.col[as.factor(expl.coll$myco.statII)])
anova(lm(HAa~coll.diff*myco.statII, data=expl.coll))
m3<-aov(HAa ~ coll.diff*myco.statII, data=expl.coll)

boxplot(HAa ~ myco.statII, data=expl.coll,col=spec.col)
TukeyHSD(m3, "myco.statII")
abline(h=0, col="gray")

## Does cons location explain PSF? YES!expl.cons
spec.col=c("green", "darkgreen", "yellow", "orange")
plot(HAa ~ sc.N.x ,data=expl.cons, pch=16,col=spec.col[as.factor(expl.cons$myco.statII)])
anova(lm(HAa ~ sc.N.x*myco.statII, data=expl.cons))
abline(lm(HAa ~ sc.N.x, data=expl.cons))
abline(h=0, col="gray")

## Does coll location explain PSF? expl.coll
spec.col[as.factor(expl.coll$growthForm.x)]
plot(HAa ~ sc.SRL.x, data=expl.coll  , pch=16,col=spec.col[as.factor(expl.coll$growthForm.x)])
anova(lm(HAa~ sc.SRL.x*myco.statII, data=expl.coll))
m3<-aov(HAa ~ sc.SRL.x*myco.statII, data=expl.coll)

boxplot(HAa ~ myco.statII, data=expl.coll,col=spec.col)
TukeyHSD(m3, "myco.statII")
abline(h=0, col="gray")
names(expl.cons)


## fancier test accounting for Species pair
names(try3)
## HA effects explained by cons diff
m1<-lmer(HAa ~ myco.statII+cons.diff *coll.diff+ (1|Family.A), data = try3) 
m2<-lmer(HAa ~ myco.statII+cons.diff +coll.diff+ (1|Family.A), data = try3)
m2b<-lmer(HAa ~ cons.diff +coll.diff+ (1|Family.A), data = try3)
anova(m1,m2)
anova(m2b,m2)

m3<-lmer(HAa ~ myco.statII+ cons.diff + (1|Family.A), data=try3) 
anova(m3,m2)
m4<-lmer(HAa ~ myco.statII+ coll.diff + (1|Family.A), data=try3) 
anova(m4,m2)

### keep this
table(try3$myco.stat)## SRL matters in AMF species 
table(try3$form.stat)## SRL becomes less important when trees or herbs are out, graminoids out= no difference
table(try3$Nfix.stat)## cannot select these

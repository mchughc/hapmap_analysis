
############
## do non-parametric sign test analysis

mxl <- get(load("hapmap_mxlOnly_estimates.RData"))
# merge in sex info
pheno <- read.csv("../HapMap3.csv",header=TRUE,as.is=TRUE)
mxl <- merge(mxl,pheno[,c("IID","sex")],by.x="V2",by.y="IID",all.x=TRUE)
sum(is.na(mxl$sex)) # 6; hmm -- looked them up online by hand
mxl$sex[is.element(mxl$V2,c("NA19672","NA19737","NA19740","NA19764"))] <- 2
mxl$sex[is.element(mxl$V2,c("NA19741","NA19792"))] <- 1
table(mxl$sex,exclude=NULL) # great
#   1    2 <NA> 
#  38   48    0 

## YRI
XAutoDiff_YRI <- mxl$chrAll.yri[mxl$unrelated]-mxl$chrX.yri[mxl$unrelated]
mean(XAutoDiff_YRI) # 0.01032438

pdf("../ancestry_differences/plots/hist_autoXDiff_Afr.pdf")
h <- hist(XAutoDiff_YRI,nclass=15,main="Difference in African Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=0.0103"),col="red",lty=2,lwd=1.5)
dev.off()


diff_ind_YRI <- XAutoDiff_YRI>0

var(mxl$chrX.yri[mxl$unrelated]) # 0.003685088
var(mxl$chrAll.yri[mxl$unrelated]) # 0.0003700711

sum(diff_ind_YRI)/length(diff_ind_YRI) # 0.6603774
# under the null, this proportion should be 1/2
m <- sum(diff_ind_YRI)
pbinom(m,53,.5) # 0.9936698
1-pbinom(m,53,.5)+dbinom(m,53,.5) # 0.01350416
2*(1-pbinom(m,53,.5)+dbinom(m,53,.5)) # 0.02700832

## stratify by sex
XAutoDiff_YRI_fem <- mxl$chrAll.yri[mxl$unrelated&mxl$sex==2]-
  mxl$chrX.yri[mxl$unrelated&mxl$sex==2]
diff_ind_YRI_fem <- XAutoDiff_YRI_fem>0

sum(diff_ind_YRI_fem)/length(diff_ind_YRI_fem) # 0.6296296
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_fem)) # 17
pbinom(m,27,.5) # 0.9389609
1-pbinom(m,27,.5)+dbinom(m,27,.5) # 0.1238943
2*(1-pbinom(m,27,.5)+dbinom(m,27,.5)) # 0.2477886

XAutoDiff_YRI_mal <- mxl$chrAll.yri[mxl$unrelated&mxl$sex==1]-
  mxl$chrX.yri[mxl$unrelated&mxl$sex==1]
diff_ind_YRI_mal <- XAutoDiff_YRI_mal>0

sum(diff_ind_YRI_mal)/length(diff_ind_YRI_mal) # 0.6923077
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_mal)) # 18
pbinom(m,26,.5) # 0.9855204
1-pbinom(m,26,.5)+dbinom(m,26,.5) # 0.03775935
2*(1-pbinom(m,26,.5)+dbinom(m,26,.5)) # 0.0755187



## CEU
XAutoDiff_CEU <- mxl$chrAll.ceu[mxl$unrelated]-mxl$chrX.ceu[mxl$unrelated]
mean(XAutoDiff_CEU) # 0.1229255

pdf("../ancestry_differences/plots/hist_autoXDiff_Euro.pdf")
h <- hist(XAutoDiff_CEU,nclass=15,main="Difference in European Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=0.1229"),col="red",lwd=1.5,lty=2)
dev.off()



diff_ind_CEU <- XAutoDiff_CEU>0

var(mxl$chrX.ceu[mxl$unrelated]) # 0.06006536
var(mxl$chrAll.ceu[mxl$unrelated]) # 0.02417326

sum(diff_ind_CEU)/length(diff_ind_CEU) # 0.7735849
m <- sum(diff_ind_CEU)
pbinom(m,53,.5) # 0.9999888
1-pbinom(m,53,.5)+dbinom(m,53,.5) # 4.085667e-05
2*(1-pbinom(m,53,.5)+dbinom(m,53,.5)) # 8.171335e-05

## stratify by sex
XAutoDiff_CEU_fem <- mxl$chrAll.ceu[mxl$unrelated&mxl$sex==2]-
  mxl$chrX.ceu[mxl$unrelated&mxl$sex==2]
diff_ind_CEU_fem <- XAutoDiff_CEU_fem>0

sum(diff_ind_CEU_fem)/length(diff_ind_CEU_fem) # 0.7777778
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_fem)) # 21
pbinom(m,27,.5) # 0.9992431
1-pbinom(m,27,.5)+dbinom(m,27,.5) # 0.002962306
2*(1-pbinom(m,27,.5)+dbinom(m,27,.5)) # 0.005924612

XAutoDiff_CEU_mal <- mxl$chrAll.ceu[mxl$unrelated&mxl$sex==1]-
  mxl$chrX.ceu[mxl$unrelated&mxl$sex==1]
diff_ind_CEU_mal <- XAutoDiff_CEU_mal>0

sum(diff_ind_CEU_mal)/length(diff_ind_CEU_mal) # 0.7692308
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_mal)) # 20
pbinom(m,26,.5) # 0.998753
1-pbinom(m,26,.5)+dbinom(m,26,.5) # 0.004677653
2*(1-pbinom(m,26,.5)+dbinom(m,26,.5)) # 0.009355307


## NAm
XAutoDiff_NAm <- mxl$chrAll.hgdp[mxl$unrelated]-mxl$chrX.hgdp[mxl$unrelated]
mean(XAutoDiff_NAm) # -0.1332499

pdf("../ancestry_differences/plots/hist_autoXDiff_NativeAm.pdf")
h <- hist(XAutoDiff_NAm,nclass=15,main="Difference in Native American Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
legend("topleft",c("mean=-0.1332"),col="red",lty=2,lwd=1.5)
dev.off()



diff_ind_NAm <- XAutoDiff_NAm>0

var(mxl$chrX.hgdp[mxl$unrelated]) # 0.0615204
var(mxl$chrAll.hgdp[mxl$unrelated]) # 0.02487759

sum(diff_ind_NAm)/length(diff_ind_NAm) # 0.2075472
# get a pvalue
m <- sum(diff_ind_NAm)
pbinom(m,53,.5) # 1.12378e-05
1-pbinom(m,53,.5)+dbinom(m,53,.5) # 0.9999972
2*pbinom(m,53,.5) # 2.247559e-05

## stratify by sex
XAutoDiff_fem <- mxl$chrAll.hgdp[mxl$unrelated&mxl$sex==2]-
  mxl$chrX.hgdp[mxl$unrelated&mxl$sex==2]
diff_ind_fem <- XAutoDiff_fem>0

sum(diff_ind_fem)/length(diff_ind_fem) # 0.1851852
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_fem)) # 5
pbinom(m,27,.5) # 0.0007568598
1-pbinom(m,27,.5)+dbinom(m,27,.5) # 0.9998446
2*pbinom(m,27,.5) # 0.00151372

XAutoDiff_mal <- mxl$chrAll.hgdp[mxl$unrelated&mxl$sex==1]-
  mxl$chrX.hgdp[mxl$unrelated&mxl$sex==1]
diff_ind_mal <- XAutoDiff_mal>0

sum(diff_ind_mal)/length(diff_ind_mal) # 0.2307692
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_mal)) # 6
pbinom(m,26,.5) # 0.004677653
1-pbinom(m,26,.5)+dbinom(m,26,.5) # 0.998753
2*pbinom(m,26,.5) # 0.009355307



#### plot of all histograms together

pdf("../ancestry_differences/plots/hist_autoXDiff_allSubpops.pdf",width=14)
par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,3))
#par(mfrow=c(2,2))


asw <- get(load("hapmap_aswOnly_estimates.RData"))
XAutoDiff_YRI <- asw$chrAll.yri[asw$unrelated]-asw$chrX.yri[asw$unrelated]
XAutoDiff_CEU <- asw$chrAll.ceu[asw$unrelated]-asw$chrX.ceu[asw$unrelated]
XAutoDiff_NAm <- asw$chrAll.hgdp[asw$unrelated]-asw$chrX.hgdp[asw$unrelated]

h <- hist(XAutoDiff_NAm,nclass=15,main="",#main="Difference in Native American Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry",cex.lab=1.5)
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
mtext("A", side=3, line=0.75,adj=0,cex=1.3)

h <- hist(XAutoDiff_CEU,nclass=15,main="",#main="Difference in European Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry",cex.lab=1.5)
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_YRI,nclass=15,main="",#main="Difference in African Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry",cex.lab=1.5)
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,lwd=1.5,col="red")

#plot(1:10,type="n",axes=F,xlab="",ylab="")
#legend("topright",c("Native Am mean difference=6.6618e-04","Euro mean difference=0.0301","African mean difference=-0.0307","x=0"),
#       lty=c(rep(2,3),1),lwd=1.5,col=c(rep("red",3),"gray"),cex=0.8)


XAutoDiff_YRI <- mxl$chrAll.yri[mxl$unrelated]-mxl$chrX.yri[mxl$unrelated]
XAutoDiff_CEU <- mxl$chrAll.ceu[mxl$unrelated]-mxl$chrX.ceu[mxl$unrelated]
XAutoDiff_NAm <- mxl$chrAll.hgdp[mxl$unrelated]-mxl$chrX.hgdp[mxl$unrelated]

h <- hist(XAutoDiff_NAm,nclass=15,main="",#main="Difference in Native American Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry",cex.lab=1.5)
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
mtext("B", side=3, line=0.75,adj=0,cex=1.3)

h <- hist(XAutoDiff_CEU,nclass=15,main="",#main="Difference in European Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry",cex.lab=1.5)
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_YRI,nclass=15,main="",#main="Difference in African Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry",cex.lab=1.5)
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,lwd=1.5,col="red")

#plot(1:10,type="n",axes=F,xlab="",ylab="")
#legend("topleft",c("Native Am mean difference=-0.1332","Euro mean difference=0.1229","African mean difference=0.0103","x=0"),
#       lty=c(rep(2,3),1),lwd=1.5,col=c(rep("red",3),"gray"),cex=0.8)

dev.off()





rm(list=ls())

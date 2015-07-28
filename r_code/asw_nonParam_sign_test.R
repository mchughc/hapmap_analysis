
############
## do non-parametric sign test analysis

asw <- get(load("hapmap_aswOnly_estimates.RData"))
# merge in sex info
pheno <- read.csv("../HapMap3.csv",header=TRUE,as.is=TRUE)
asw <- merge(asw,pheno[,c("IID","sex")],by.x="V2",by.y="IID",all.x=TRUE)
sum(is.na(asw$sex)) # 4; hmm -- looked them up online by hand
asw$sex[is.element(asw$V2,c("NA19984","NA20298","NA20351"))] <- 1
asw$sex[is.element(asw$V2,"NA20412")] <- 2
table(asw$sex,exclude=NULL) # great
#   1    2 <NA> 
#  41   46    0 

## YRI
XAutoDiff_YRI <- asw$chrAll.yri[asw$unrelated]-asw$chrX.yri[asw$unrelated]
mean(XAutoDiff_YRI) # -0.03071591

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_Afr.pdf")
h <- hist(XAutoDiff_YRI,nclass=15,main="Difference in African Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=-0.0307"),col="red",lty=2,lwd=1.5)
dev.off()


diff_ind_YRI <- XAutoDiff_YRI>0

var(asw$chrX.yri[asw$unrelated]) # 0.02337575
var(asw$chrAll.yri[asw$unrelated]) # 0.006532156

sum(diff_ind_YRI)/length(diff_ind_YRI) # 0.3777778
# under the null, this proportion should be 1/2
m <- sum(diff_ind_YRI) # 17
pbinom(m,45,.5) # 0.06757823
2*pbinom(m,45,.5) # 0.1351565

## stratify by sex
XAutoDiff_YRI_fem <- asw$chrAll.yri[asw$unrelated&asw$sex==2]-
  asw$chrX.yri[asw$unrelated&asw$sex==2]
diff_ind_YRI_fem <- XAutoDiff_YRI_fem>0

sum(diff_ind_YRI_fem)/length(diff_ind_YRI_fem) # 0.44
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_fem)) # 11
pbinom(m,25,.5) # 0.345019
2*pbinom(m,25,.5) # 0.690038
1-pbinom(m,45,.5)+dbinom(m,45,.5) # 0.03622713
2*(1-pbinom(m,45,.5)+dbinom(m,45,.5)) # 0.07245426

XAutoDiff_YRI_mal <- asw$chrAll.yri[asw$unrelated&asw$sex==1]-
  asw$chrX.yri[asw$unrelated&asw$sex==1]
diff_ind_YRI_mal <- XAutoDiff_YRI_mal>0

sum(diff_ind_YRI_mal)/length(diff_ind_YRI_mal) # 0.3
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_mal)) # 6
pbinom(m,20,.5) # 0.05765915
2*pbinom(m,20,.5) # 0.1153183



## CEU
XAutoDiff_CEU <- asw$chrAll.ceu[asw$unrelated]-asw$chrX.ceu[asw$unrelated]
mean(XAutoDiff_CEU) # 0.03005403

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_Euro.pdf")
h <- hist(XAutoDiff_CEU,nclass=15,main="Difference in European Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=0.0301"),col="red",lwd=1.5,lty=2)
dev.off()



diff_ind_CEU <- XAutoDiff_CEU>0

var(asw$chrX.ceu[asw$unrelated]) # 0.02060563
var(asw$chrAll.ceu[asw$unrelated]) # 0.005607704

sum(diff_ind_CEU)/length(diff_ind_CEU) # 0.6444444
m <- sum(diff_ind_CEU) # 29
pbinom(m,45,.5) # 0.9821511
1-pbinom(m,45,.5)+dbinom(m,45,.5) # 0.03622713
2*(1-pbinom(m,45,.5)+dbinom(m,45,.5)) # 0.07245426

## stratify by sex
XAutoDiff_CEU_fem <- asw$chrAll.ceu[asw$unrelated&asw$sex==2]-
  asw$chrX.ceu[asw$unrelated&asw$sex==2]
diff_ind_CEU_fem <- XAutoDiff_CEU_fem>0

sum(diff_ind_CEU_fem)/length(diff_ind_CEU_fem) # 0.56
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_fem)) # 14
pbinom(m,25,.5) # 0.7878219
1-pbinom(m,25,.5)+dbinom(m,25,.5) # 0.345019
2*(1-pbinom(m,25,.5)+dbinom(m,25,.5)) # 0.690038

XAutoDiff_CEU_mal <- asw$chrAll.ceu[asw$unrelated&asw$sex==1]-
  asw$chrX.ceu[asw$unrelated&asw$sex==1]
diff_ind_CEU_mal <- XAutoDiff_CEU_mal>0

sum(diff_ind_CEU_mal)/length(diff_ind_CEU_mal) # 0.75
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_mal)) # 15
pbinom(m,20,.5) # 0.994091
1-pbinom(m,20,.5)+dbinom(m,20,.5) # 0.02069473
2*(1-pbinom(m,20,.5)+dbinom(m,20,.5)) # 0.04138947


## NAm
XAutoDiff_NAm <- asw$chrAll.hgdp[asw$unrelated]-asw$chrX.hgdp[asw$unrelated]
mean(XAutoDiff_NAm) # 0.0006617879

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_NativeAm.pdf")
h <- hist(XAutoDiff_NAm,nclass=15,main="Difference in Native American Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
legend("topleft",c("mean=6.62e-04"),col="red",lty=2,lwd=1.5)
dev.off()



diff_ind_NAm <- XAutoDiff_NAm>0

var(asw$chrX.hgdp[asw$unrelated]) # 0.00127514
var(asw$chrAll.hgdp[asw$unrelated]) # 0.001000164

sum(diff_ind_NAm)/length(diff_ind_NAm) # 0.6666667
# get a pvalue
m <- sum(diff_ind_NAm) # 30
pbinom(m,45,.5) # 0.9919528
1-pbinom(m,45,.5)+dbinom(m,45,.5) # 0.0178489
2*(1-pbinom(m,45,.5)+dbinom(m,45,.5)) # 0.0356978

## stratify by sex
XAutoDiff_fem <- asw$chrAll.hgdp[asw$unrelated&asw$sex==2]-
  asw$chrX.hgdp[asw$unrelated&asw$sex==2]
diff_ind_fem <- XAutoDiff_fem>0

sum(diff_ind_fem)/length(diff_ind_fem) # 0.76
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_fem)) # 19
pbinom(m,25,.5) # 0.9979613
1-pbinom(m,25,.5)+dbinom(m,25,.5) # 0.007316649
2*(1-pbinom(m,25,.5)+dbinom(m,25,.5)) # 0.0146333

XAutoDiff_mal <- asw$chrAll.hgdp[asw$unrelated&asw$sex==1]-
  asw$chrX.hgdp[asw$unrelated&asw$sex==1]
diff_ind_mal <- XAutoDiff_mal>0

sum(diff_ind_mal)/length(diff_ind_mal) # 0.55
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_mal)) # 11
pbinom(m,20,.5) # 0.7482777
1-pbinom(m,20,.5)+dbinom(m,20,.5) # 0.4119015
2*(1-pbinom(m,20,.5)+dbinom(m,20,.5)) # 0.8238029


#### plot of all histograms together

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_allSubpops.pdf")
par(mfrow=c(2,2))
h <- hist(XAutoDiff_NAm,nclass=15,main="Difference in Native American Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_CEU,nclass=15,main="Difference in European Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_YRI,nclass=15,main="Difference in African Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,lwd=1.5,col="red")

plot(1:10,type="n",axes=F,xlab="",ylab="")
legend("topright",c("Native Am mean difference=6.6618e-04","Euro mean difference=0.0301","African mean difference=-0.0307","x=0"),
       lty=c(rep(2,3),1),lwd=1.5,col=c(rep("red",3),"gray"),cex=0.8)
dev.off()





rm(list=ls())

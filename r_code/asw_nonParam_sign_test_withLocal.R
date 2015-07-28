
### non param sign test in asw samples using local ancestry estimates


############
## do non-parametric sign test analysis

asw <- get(load("asw_estimates_withLocal.RData"))
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
XAutoDiff_YRI <- asw$chrAll.local.YRI[asw$unrelated]-asw$chrX.local.YRI[asw$unrelated]
mean(XAutoDiff_YRI) # -0.03314869

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_Afr_localAncest.pdf")
h <- hist(XAutoDiff_YRI,nclass=15,main="Difference in Avg Local African Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=-0.0331"),col="red",lty=2,lwd=1.5)
dev.off()


diff_ind_YRI <- XAutoDiff_YRI>0

var(asw$chrX.local.YRI[asw$unrelated]) # 0.01951959
var(asw$chrAll.local.YRI[asw$unrelated]) # 0.006803858

sum(diff_ind_YRI)/length(diff_ind_YRI) # 0.4
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI)) # 18
pbinom(m,45,.5) # 0.1163466
2*pbinom(m,45,.5) # 0.2326932

## stratify by sex
XAutoDiff_YRI_fem <- asw$chrAll.local.YRI[asw$unrelated&asw$sex==2]-
  asw$chrX.local.YRI[asw$unrelated&asw$sex==2]
diff_ind_YRI_fem <- XAutoDiff_YRI_fem>0

sum(diff_ind_YRI_fem)/length(diff_ind_YRI_fem) # 0.48
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_fem)) # 12
pbinom(m,25,.5) # 0.5
2*pbinom(m,25,.5) # 1
1-pbinom(m,45,.5)+dbinom(m,45,.5) # 0.999588
2*(1-pbinom(m,45,.5)+dbinom(m,45,.5)) # 1.999176

XAutoDiff_YRI_mal <- asw$chrAll.local.YRI[asw$unrelated&asw$sex==1]-
  asw$chrX.local.YRI[asw$unrelated&asw$sex==1]
diff_ind_YRI_mal <- XAutoDiff_YRI_mal>0

sum(diff_ind_YRI_mal)/length(diff_ind_YRI_mal) # 0.3
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_mal)) # 6
pbinom(m,20,.5) # 0.05765915
2*pbinom(m,20,.5) # 0.1153183



## CEU
XAutoDiff_CEU <- asw$chrAll.local.CEU[asw$unrelated]-asw$chrX.local.CEU[asw$unrelated]
mean(XAutoDiff_CEU) # 0.03340967

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_Euro_localAncest.pdf")
h <- hist(XAutoDiff_CEU,nclass=15,main="Difference in Avg Local European Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=0.0334"),col="red",lwd=1.5,lty=2)
dev.off()



diff_ind_CEU <- XAutoDiff_CEU>0

var(asw$chrX.local.CEU[asw$unrelated]) # 0.01915766
var(asw$chrAll.local.CEU[asw$unrelated]) # 0.005990838

sum(diff_ind_CEU)/length(diff_ind_CEU) # 0.6222222
(m <- sum(diff_ind_CEU)) # 28
pbinom(m,45,.5) # 0.9637729
1-pbinom(m,45,.5)+dbinom(m,45,.5) # 0.06757823
2*(1-pbinom(m,45,.5)+dbinom(m,45,.5)) # 0.1351565

## stratify by sex
XAutoDiff_CEU_fem <- asw$chrAll.local.CEU[asw$unrelated&asw$sex==2]-
  asw$chrX.local.CEU[asw$unrelated&asw$sex==2]
diff_ind_CEU_fem <- XAutoDiff_CEU_fem>0

sum(diff_ind_CEU_fem)/length(diff_ind_CEU_fem) # 0.56
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_fem)) # 14
pbinom(m,25,.5) # 0.7878219
1-pbinom(m,25,.5)+dbinom(m,25,.5) # 0.345019
2*(1-pbinom(m,25,.5)+dbinom(m,25,.5)) # 0.690038

XAutoDiff_CEU_mal <- asw$chrAll.local.CEU[asw$unrelated&asw$sex==1]-
  asw$chrX.local.CEU[asw$unrelated&asw$sex==1]
diff_ind_CEU_mal <- XAutoDiff_CEU_mal>0

sum(diff_ind_CEU_mal)/length(diff_ind_CEU_mal) # 0.7
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_mal)) # 14
pbinom(m,20,.5) # 0.9793053
1-pbinom(m,20,.5)+dbinom(m,20,.5) # 0.05765915
2*(1-pbinom(m,20,.5)+dbinom(m,20,.5)) # 0.1153183


## NAm
XAutoDiff_NAm <- asw$chrAll.local.NAM[asw$unrelated]-asw$chrX.local.NAM[asw$unrelated]
mean(XAutoDiff_NAm) # -0.0002609808

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_NativeAm_localAncest.pdf")
h <- hist(XAutoDiff_NAm,nclass=15,main="Difference in Avg Local Native American Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
legend("topleft",c("mean=-2.61e-04"),col="red",lty=2,lwd=1.5)
dev.off()



diff_ind_NAm <- XAutoDiff_NAm>0

var(asw$chrX.local.NAM[asw$unrelated]) # 0.001340183
var(asw$chrAll.local.NAM[asw$unrelated]) # 0.001024059

sum(diff_ind_NAm)/length(diff_ind_NAm) # 0.7555556
# get a pvalue
(m <- sum(diff_ind_NAm)) # 34
pbinom(m,45,.5) # 0.9998765
1-pbinom(m,45,.5)+dbinom(m,45,.5) # 0.0004120412
2*(1-pbinom(m,45,.5)+dbinom(m,45,.5)) # 0.0008240824

## stratify by sex
XAutoDiff_fem <- asw$chrAll.local.NAM[asw$unrelated&asw$sex==2]-
  asw$chrX.local.NAM[asw$unrelated&asw$sex==2]
diff_ind_fem <- XAutoDiff_fem>0

sum(diff_ind_fem)/length(diff_ind_fem) # 0.68
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_fem)) # 17
pbinom(m,25,.5) # 0.9783574
1-pbinom(m,25,.5)+dbinom(m,25,.5) # 0.05387607
2*(1-pbinom(m,25,.5)+dbinom(m,25,.5)) # 0.1077521

XAutoDiff_mal <- asw$chrAll.local.NAM[asw$unrelated&asw$sex==1]-
  asw$chrX.local.NAM[asw$unrelated&asw$sex==1]
diff_ind_mal <- XAutoDiff_mal>0

sum(diff_ind_mal)/length(diff_ind_mal) # 0.85
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_mal)) # 17
pbinom(m,20,.5) # 0.9997988
1-pbinom(m,20,.5)+dbinom(m,20,.5) # 0.001288414
2*(1-pbinom(m,20,.5)+dbinom(m,20,.5)) # 0.002576828


#### plot of all histograms together

pdf("../ancestry_differences/plots/asw_hist_autoXDiff_allSubpops_localAncest.pdf")
par(mfrow=c(2,2))
h <- hist(XAutoDiff_NAm,nclass=15,main="Difference in Avg Local Native American Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_CEU,nclass=15,main="Difference in Avg Local European Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_YRI,nclass=15,main="Difference in Avg Local African Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,lwd=1.5,col="red")

plot(1:10,type="n",axes=F,xlab="",ylab="")
legend("topright",c("Native Am mean difference=-2.61e-04","Euro mean difference=0.0334","African mean difference=-0.0331","x=0"),
       lty=c(rep(2,3),1),lwd=1.5,col=c(rep("red",3),"gray"),cex=0.8)
dev.off()


#### get corr estimates by chr comparing the local and ADMIX estimates
localEst <- asw[,79:147]
admixEst <- asw[,4:72]
corValues <- rep(NA,69)

for(i in 1:ncol(localEst)){
  print(colnames(localEst)[i])
  print(cor(localEst[,i],admixEst[,i]))
  corValues[i] <- cor(localEst[,i],admixEst[,i])
  names(corValues[i]) <- colnames(localEst)[i]
}

# autosomal correlation
localEst <- asw[,146:148]
admixEst <- asw[,73:75]
corValues <- rep(NA,3)

for(i in 1:ncol(localEst)){
  corValues[i] <- cor(localEst[,i],admixEst[,i])
}
corValues

# x chr correlation
localEst <- asw[,143:145]
admixEst <- asw[,70:72]
corValues <- rep(NA,3)

for(i in 1:ncol(localEst)){
  corValues[i] <- cor(localEst[,i],admixEst[,i])
}
corValues





rm(list=ls())



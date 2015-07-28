
# non param sign test on local ancestry results in MXL samples




############
## do non-parametric sign test analysis

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
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
XAutoDiff_YRI <- mxl$chrAll.local.YRI[mxl$unrelated]-mxl$chrX.local.YRI[mxl$unrelated]
mean(XAutoDiff_YRI) # 0.003259286

pdf("../ancestry_differences/plots/hist_autoXDiff_Afr_localAncest.pdf")
h <- hist(XAutoDiff_YRI,nclass=15,main="Difference in Avg Local African Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_YRI),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=0.0033"),col="red",lty=2,lwd=1.5)
dev.off()


diff_ind_YRI <- XAutoDiff_YRI>0

var(mxl$chrX.local.YRI[mxl$unrelated]) # 0.003554825
var(mxl$chrAll.local.YRI[mxl$unrelated]) # 0.0003723114

sum(diff_ind_YRI)/length(diff_ind_YRI) # 0.6415094
# under the null, this proportion should be 1/2
m <- sum(diff_ind_YRI)
pbinom(m,53,.5) # 0.9864958
1-pbinom(m,53,.5)+dbinom(m,53,.5) # 0.02671941
2*(1-pbinom(m,53,.5)+dbinom(m,53,.5)) # 0.05343881

## stratify by sex
XAutoDiff_YRI_fem <- mxl$chrAll.local.YRI[mxl$unrelated&mxl$sex==2]-
  mxl$chrX.local.YRI[mxl$unrelated&mxl$sex==2]
diff_ind_YRI_fem <- XAutoDiff_YRI_fem>0

sum(diff_ind_YRI_fem)/length(diff_ind_YRI_fem) # 0.5925926
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_fem)) # 16
pbinom(m,27,.5) # 0.8761057
1-pbinom(m,27,.5)+dbinom(m,27,.5) # 0.2210342
2*(1-pbinom(m,27,.5)+dbinom(m,27,.5)) # 0.4420683

XAutoDiff_YRI_mal <- mxl$chrAll.local.YRI[mxl$unrelated&mxl$sex==1]-
  mxl$chrX.local.YRI[mxl$unrelated&mxl$sex==1]
diff_ind_YRI_mal <- XAutoDiff_YRI_mal>0

sum(diff_ind_YRI_mal)/length(diff_ind_YRI_mal) # 0.6923077
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_YRI_mal)) # 18
pbinom(m,26,.5) # 0.9855204
1-pbinom(m,26,.5)+dbinom(m,26,.5) # 0.03775935
2*(1-pbinom(m,26,.5)+dbinom(m,26,.5)) # 0.0755187



## CEU
XAutoDiff_CEU <- mxl$chrAll.local.CEU[mxl$unrelated]-mxl$chrX.local.CEU[mxl$unrelated]
mean(XAutoDiff_CEU) # 0.1290034

pdf("../ancestry_differences/plots/hist_autoXDiff_Euro_localAncest.pdf")
h <- hist(XAutoDiff_CEU,nclass=15,main="Difference in Avg Local European Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_CEU),lty=2,col="red",lwd=1.5)
legend("topleft",c("mean=0.1290"),col="red",lwd=1.5,lty=2)
dev.off()



diff_ind_CEU <- XAutoDiff_CEU>0

var(mxl$chrX.local.CEU[mxl$unrelated]) # 0.05794262
var(mxl$chrAll.local.CEU[mxl$unrelated]) # 0.02374341

sum(diff_ind_CEU)/length(diff_ind_CEU) # 0.8113208
(m <- sum(diff_ind_CEU)) # 43
pbinom(m,53,.5) # 0.9999994
1-pbinom(m,53,.5)+dbinom(m,53,.5) # 2.77526e-06
2*(1-pbinom(m,53,.5)+dbinom(m,53,.5)) # 5.550521e-06

## stratify by sex
XAutoDiff_CEU_fem <- mxl$chrAll.local.CEU[mxl$unrelated&mxl$sex==2]-
  mxl$chrX.local.CEU[mxl$unrelated&mxl$sex==2]
diff_ind_CEU_fem <- XAutoDiff_CEU_fem>0

sum(diff_ind_CEU_fem)/length(diff_ind_CEU_fem) # 0.8148148
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_fem)) # 22
pbinom(m,27,.5) # 0.9998446
1-pbinom(m,27,.5)+dbinom(m,27,.5) # 0.0007568598
2*(1-pbinom(m,27,.5)+dbinom(m,27,.5)) # 0.00151372

XAutoDiff_CEU_mal <- mxl$chrAll.local.CEU[mxl$unrelated&mxl$sex==1]-
  mxl$chrX.local.CEU[mxl$unrelated&mxl$sex==1]
diff_ind_CEU_mal <- XAutoDiff_CEU_mal>0

sum(diff_ind_CEU_mal)/length(diff_ind_CEU_mal) # 0.8076923
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_CEU_mal)) # 21
pbinom(m,26,.5) # 0.9997332
1-pbinom(m,26,.5)+dbinom(m,26,.5) # 0.001246959
2*(1-pbinom(m,26,.5)+dbinom(m,26,.5)) # 0.002493918


## NAm
XAutoDiff_NAm <- mxl$chrAll.local.NAM[mxl$unrelated]-mxl$chrX.local.NAM[mxl$unrelated]
mean(XAutoDiff_NAm) # -0.1322627

pdf("../ancestry_differences/plots/hist_autoXDiff_NativeAm_localAncest.pdf")
h <- hist(XAutoDiff_NAm,nclass=15,main="Difference in Avg Local Native American Ancestry\nBetween Autosomes and X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry")
abline(v=0,col="gray")
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], lwd = 2)
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
legend("topleft",c("mean=-0.1323"),col="red",lty=2,lwd=1.5)
dev.off()



diff_ind_NAm <- XAutoDiff_NAm>0

var(mxl$chrX.local.NAM[mxl$unrelated]) # 0.05711534
var(mxl$chrAll.local.NAM[mxl$unrelated]) # 0.02357721

sum(diff_ind_NAm)/length(diff_ind_NAm) # 0.1886792
# get a pvalue
(m <- sum(diff_ind_NAm)) # 10
pbinom(m,53,.5) # 2.77526e-06
1-pbinom(m,53,.5)+dbinom(m,53,.5) # 0.9999994
2*pbinom(m,53,.5) # 5.550521e-06

## stratify by sex
XAutoDiff_fem <- mxl$chrAll.local.NAM[mxl$unrelated&mxl$sex==2]-
  mxl$chrX.local.NAM[mxl$unrelated&mxl$sex==2]
diff_ind_fem <- XAutoDiff_fem>0

sum(diff_ind_fem)/length(diff_ind_fem) # 0.1851852
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_fem)) # 5
pbinom(m,27,.5) # 0.0007568598
1-pbinom(m,27,.5)+dbinom(m,27,.5) # 0.9998446
2*pbinom(m,27,.5) # 0.00151372

XAutoDiff_mal <- mxl$chrAll.local.NAM[mxl$unrelated&mxl$sex==1]-
  mxl$chrX.local.NAM[mxl$unrelated&mxl$sex==1]
diff_ind_mal <- XAutoDiff_mal>0

sum(diff_ind_mal)/length(diff_ind_mal) # 0.1923077
# under the null, this proportion should be 1/2
(m <- sum(diff_ind_mal)) # 5
pbinom(m,26,.5) #  0.001246959
1-pbinom(m,26,.5)+dbinom(m,26,.5) # 0.9997332
2*pbinom(m,26,.5) # 0.002493918



#### plot of all histograms together

pdf("../ancestry_differences/plots/hist_autoXDiff_allSubpops_localAncest.pdf",width=14)
par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,3))
#par(mfrow=c(2,2))

asw <- get(load("asw_estimates_withLocal.RData"))
XAutoDiff_YRI <- asw$chrAll.local.YRI[asw$unrelated]-asw$chrX.local.YRI[asw$unrelated]
XAutoDiff_CEU <- asw$chrAll.local.CEU[asw$unrelated]-asw$chrX.local.CEU[asw$unrelated]
XAutoDiff_NAm <- asw$chrAll.local.NAM[asw$unrelated]-asw$chrX.local.NAM[asw$unrelated]

h <- hist(XAutoDiff_NAm,nclass=15,main="",#main="Difference in Avg Local Native American Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry",cex.lab=1.5)
abline(v=0,lwd=2)
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], col="gray")
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
mtext("A", side=3, line=0.75,adj=0,cex=1.3)

h <- hist(XAutoDiff_CEU,nclass=15,main="",#main="Difference in Avg Local European Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry",cex.lab=1.5)
abline(v=0,lwd=2)
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], col="gray")
abline(v=mean(XAutoDiff_CEU),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_YRI,nclass=15,main="",#main="Difference in Avg Local African Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry",cex.lab=1.5)
abline(v=0,lwd=2)
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], col="gray")
abline(v=mean(XAutoDiff_YRI),lty=2,lwd=1.5,col="red")

#plot(1:10,type="n",axes=F,xlab="",ylab="")
#legend("topright",c("Native Am mean difference=6.6618e-04","Euro mean difference=0.0301","African mean difference=-0.0307","x=0"),
#       lty=c(rep(2,3),1),lwd=1.5,col=c(rep("red",3),"gray"),cex=0.8)



XAutoDiff_YRI <- mxl$chrAll.local.YRI[mxl$unrelated]-mxl$chrX.local.YRI[mxl$unrelated]
XAutoDiff_CEU <- mxl$chrAll.local.CEU[mxl$unrelated]-mxl$chrX.local.CEU[mxl$unrelated]
XAutoDiff_NAm <- mxl$chrAll.local.NAM[mxl$unrelated]-mxl$chrX.local.NAM[mxl$unrelated]

h <- hist(XAutoDiff_NAm,nclass=15,main="",#main="Difference in Avg Local Native American Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) Native American Ancestry",cex.lab=1.5)
abline(v=0,lwd=2)
d <- density(XAutoDiff_NAm)
lines(x = d$x, y = d$y * length(XAutoDiff_NAm) * diff(h$breaks)[1], col="gray")
abline(v=mean(XAutoDiff_NAm),lty=2,lwd=1.5,col="red")
mtext("B", side=3, line=0.75,adj=0,cex=1.3)

h <- hist(XAutoDiff_CEU,nclass=15,main="",#main="Difference in Avg Local European Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) European Ancestry",cex.lab=1.5)
abline(v=0,lwd=2)
d <- density(XAutoDiff_CEU)
lines(x = d$x, y = d$y * length(XAutoDiff_CEU) * diff(h$breaks)[1], col="gray")
abline(v=mean(XAutoDiff_CEU),lty=2,lwd=1.5,col="red")

h <- hist(XAutoDiff_YRI,nclass=15,main="",#main="Difference in Avg Local African Ancestry\nBetween Autosomes & X Chromosome",
          xlab="(Autosomal - X Chr) African Ancestry",cex.lab=1.5)
abline(v=0,lwd=2)
d <- density(XAutoDiff_YRI)
lines(x = d$x, y = d$y * length(XAutoDiff_YRI) * diff(h$breaks)[1], col="gray")
abline(v=mean(XAutoDiff_YRI),lty=2,lwd=1.5,col="red")

#plot(1:10,type="n",axes=F,xlab="",ylab="")
#legend("topleft",c("Native Am mean difference=-0.1332","Euro mean difference=0.1229","African mean difference=0.0103","x=0"),
#       lty=c(rep(2,3),1),lwd=1.5,col=c(rep("red",3),"gray"),cex=0.8)

dev.off()



#### get corr estimates by chr comparing the local and ADMIX estimates
localEst <- mxl[,79:147]
admixEst <- mxl[,4:72]
corValues <- rep(NA,69)

for(i in 1:ncol(localEst)){
  print(colnames(localEst)[i])
  print(cor(localEst[,i],admixEst[,i]))
  corValues[i] <- cor(localEst[,i],admixEst[,i])
  names(corValues[i]) <- colnames(localEst)[i]
}
  
# autosomal correlation
localEst <- mxl[,148:150]
admixEst <- mxl[,73:75]
corValues <- rep(NA,3)

for(i in 1:ncol(localEst)){
  corValues[i] <- cor(localEst[,i],admixEst[,i])
}

# x chr correlation
localEst <- mxl[,145:147]
admixEst <- mxl[,70:72]
corValues <- rep(NA,3)

for(i in 1:ncol(localEst)){
  corValues[i] <- cor(localEst[,i],admixEst[,i])
}


# make a plot of the table 3 pvalues, for mxl and asw independently
# 3 dots at each ancestry for all samples, just males, and just females
# xaxis=afr, euro, Nam
# yaxis=-log10(pvalue)

res <- data.frame(Ancestry=rep(c("African","European","Native American"),times=3),set=c(rep("All Samples",3),rep("Females",3),rep("Males",3)))
res <- rbind(res,res)
res$HapMap <- c(rep("MXL",9),rep("ASW",9))

# so first are MXL; all; Afr, Euro, NAm
res$pvalue <- c(0.05343881,5.550521e-06, 0.002493918, # afr, euro, nam in all MXL samp                
                0.4420683, 0.00151372, 0.00151372, # afr, euro, nam in female MXL samp
                0.0755187, 0.002493918, 0.002493918, # afr euro nam in male MXL samp
                0.2326932,0.1351565,0.0008240824, # afr euro nam in all ASW samp
                0.5,0.690038,0.1077521, # afr euro nam in female ASW samp
                0.1153183, 0.1153183, 0.002576828 # afr euro nam in male ASW samp
)

library(ggplot2)
pdf("../ancestry_differences/plots/non_param_cand_pvalues_localAncest.pdf")
ggplot(res,aes(x=Ancestry,y=-log10(pvalue),color=set,shape=HapMap)) + geom_point()
dev.off()


rm(list=ls())



## hapmap MXL results using local ancestry estimates


## load in all data
mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
dim(mxl) # 86 150


# check proportion native american --
png("mxl_chrall_chrX_hgdp_localAncest.png")
plot(mxl[,"chrAll.local.NAM"],mxl[,"chrX.local.NAM"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of Avg Local HGDP Americas Ancestry")
points(mxl[mxl$unrelated,"chrAll.local.NAM"],mxl[mxl$unrelated,"chrX.local.NAM"],col="red")
points(mxl[!mxl$unrelated,"chrAll.local.NAM"],mxl[!mxl$unrelated,"chrX.local.NAM"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
dev.off()

# make a barplot of the x chr results
png("mxl_chrx_frappe_localAncest.png")
xchr <- mxl[,c("chrX.local.CEU","chrX.local.YRI","chrX.local.NAM")]
xchr <- xchr[order(xchr$chrX.local.CEU,xchr$chrX.local.NAM,xchr$chrX.local.YRI),]
tt <- t(xchr)
colnames(tt) <- rep("",ncol(tt))
barplot(tt,col=c("blue","red","green"),xlab="Individual",main="HapMap MXL Estimated Avg Local Ancestry\nX Chromosome")
dev.off()

l <- lm(mxl[mxl$unrelated,"chrX.local.NAM"]~mxl[mxl$unrelated,"chrAll.local.NAM"])
png("mxl_chrall_chrX_hgdp_unrel_localAncest.png")
plot(mxl[,"chrAll.local.NAM"],mxl[,"chrX.local.NAM"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of Avg Local HGDP Americas Ancestry\nFor 53 Unrelated HapMap MXL")
points(mxl[mxl$unrelated,"chrAll.local.NAM"],mxl[mxl$unrelated,"chrX.local.NAM"],col="black")
#points(mxl[!mxl$unrelated,"chrAll.local.NAM"],mxl[!mxl$unrelated,"chrX.local.NAM"],col="black")
#legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
abline(summary(l)$coef[1,1],summary(l)$coef[2,1])
dev.off()

var(mxl$chrX.local.NAM[mxl$unrelated]) #  0.05711534
var(mxl$chrAll.local.NAM[mxl$unrelated]) # 0.02357721
t.test(mxl$chrX.local.NAM[mxl$unrelated],mxl$chrAll.local.NAM[mxl$unrelated])
# t = 3.3897, df = 88.681, p-value = 0.001046
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.05472841 0.20979690
# sample estimates:
#   mean of x mean of y 
# 0.5812864 0.4490237 


var(mxl$chr8.local.NAM[mxl$unrelated]) # 0.02853808
var(mxl$chr15.local.NAM[mxl$unrelated]) # 0.03354886
t.test(mxl$chr8.local.NAM[mxl$unrelated],mxl$chr15.local.NAM[mxl$unrelated])
# t = 0.2922, df = 103.327, p-value = 0.7707
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.05787563  0.07787955
# sample estimates:
#   mean of x mean of y 
# 0.4494281 0.4394262 

var(mxl$chrX.local.NAM[mxl$unrelated]) #  0.05711534
var(mxl$chr15.local.NAM[mxl$unrelated]) # 0.03354886
t.test(mxl$chrX.local.NAM[mxl$unrelated],mxl$chr15.local.NAM[mxl$unrelated])
# t = 3.4299, df = 97.418, p-value = 0.0008869
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.05977663 0.22394383
# sample estimates:
#   mean of x mean of y 
# 0.5812864 0.4394262 

##################
# do analysis comparing each chr to the pool of all other chrs
res <- data.frame(matrix(NA,nrow=23,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=81,to=147,by=3)]
colnames(sm)
for(i in 1:23){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*23)
}       
res # xchr bonf pvalue: 2.612223e-05

resEuro <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=79,to=145,by=3)]       
colnames(sm)
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*23)
}       
resEuro # xchr bonf pvalue: 2.109122e-05

resAfr <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=80,to=146,by=3)]       
colnames(sm)
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resAfr[i,] <- c(i,t$p.value,t$p.value*23)
}       
resAfr # nothing!

library(xtable)
xtable(cbind(res[,2:3],resEuro[,2:3],resAfr[,2:3]))

resAfr$bonf_pvalue[resAfr$bonf_pvalue>1]<-1
resEuro$bonf_pvalue[resEuro$bonf_pvalue>1]<-1
res$bonf_pvalue[res$bonf_pvalue>1]<-1
xtable(cbind(res[,2:3],resEuro[,2:3],resAfr[,2:3]))

###
# combine the p-values using fisher's combined probability test
(chisq_cand <- -2*sum(log(res$pvalue))) # 83.20645
(afr_chisq_cand <- -2*sum(log(resAfr$pvalue))) # 43.2722
(euro_chisq_cand <- -2*sum(log(resEuro$pvalue))) # 77.03308

1-pchisq(chisq_cand,df=2*23) # 0.0006431172
1-pchisq(afr_chisq_cand,df=2*23) # 0.5871784
1-pchisq(euro_chisq_cand,df=2*23) # 0.002793493

#png("paired_ttest_pools.png")
pdf("paired_ttest_pools_localAncest.pdf")
plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="P-values from Paired T-Tests")
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.8)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"))
dev.off()

#png("paired_ttest_bonfCorr_pools.png")
pdf("paired_ttest_bonfCorr_pools_localAncest.pdf")
plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="Bonferroni Corrected P-values from Paired T-Tests")
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.8)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"))
dev.off()


#png("paired_ttest_bonfAndNot.png")
pdf("../ancestry_differences/plots/paired_ttest_bonfAndNot_paper_localAncest.pdf")
par( mai=c(0.65, 0.65, 0.4, 0.18), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))

plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE)
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)

mtext("A", side=3, line=0.75,adj=0,cex=1.3)

plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE)
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
mtext("B", side=3, line=0.75,adj=0,cex=1.3)

dev.off()



#################
# do analysis comparing each chr to the pool of all other chrs
# excluding the X chromosome
res <- data.frame(matrix(NA,nrow=22,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=81,to=144,by=3)]
colnames(sm)
for(i in 1:22){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*22)
}       

resEuro <- data.frame(matrix(NA,nrow=22,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=79,to=142,by=3)]       
colnames(sm)
for(i in 1:22){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*22)
}       

resAfr <- data.frame(matrix(NA,nrow=22,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=80,to=143,by=3)]       
colnames(sm)
for(i in 1:22){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resAfr[i,] <- c(i,t$p.value,t$p.value*22)
}       
resAfr # nothing!

library(xtable)
resAfr$bonf_pvalue[resAfr$bonf_pvalue>1]<-1
resEuro$bonf_pvalue[resEuro$bonf_pvalue>1]<-1
res$bonf_pvalue[res$bonf_pvalue>1]<-1
xtable(cbind(res[,2:3],resEuro[,2:3],resAfr[,2:3]))

bigMax <- max(c(-log10(res$pvalue),-log10(resEuro$pvalue),-log10(resAfr$pvalue)))
#png("paired_ttest_pools_exclXchr.png")
pdf("../ancestry_differences/plots/paired_ttest_pools_exclXchr_localAncest.pdf")
plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="P-values from Paired T-Tests",ylim=c(0,bigMax))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
axis(1,at=1:22,labels=c(1:22),cex.axis=0.8)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)
dev.off()

bigMax <- max(c(-log10(res$bonf_pvalue),-log10(resEuro$bonf_pvalue),-log10(resAfr$bonf_pvalue)))
#png("paired_ttest_bonfCorr_pools_exclXchr.png")
pdf("../ancestry_differences/plots/paired_ttest_bonfCorr_pools_exclXchr_localAncest.pdf")
plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="Bonferroni Corrected P-values from Paired T-Tests",ylim=c(0,bigMax))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
axis(1,at=1:22,labels=c(1:22),cex.axis=0.8)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)
dev.off()

#png("paired_ttest_bonfAndNot_exclXchr.png")
pdf("../ancestry_differences/plots/paired_ttest_bonfAndNot_exclXchr_localAncest.pdf")
par( mai=c(0.65, 0.65, 0.4, 0.15), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
bigMax <- max(c(-log10(res$pvalue),-log10(resEuro$pvalue),-log10(resAfr$pvalue)))
plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     ylim=c(0,bigMax))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:22,labels=c(1:22),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)
bigMax <- max(c(-log10(res$bonf_pvalue),-log10(resEuro$bonf_pvalue),-log10(resAfr$bonf_pvalue)))
plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     ylim=c(0,2))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:22,labels=c(1:22),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topright",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)
mtext("B", side=3, line=0.75,adj=0,cex=1.3)
dev.off()


#### make boxplots of ancestry proportions for each ancestry subpop
pdf("../ancestry_differences/plots/boxplot_ancestryProp_localAncest.pdf",width=10)
cols <- c("chrAll.local.CEU","chrX.local.CEU","chrAll.local.YRI","chrX.local.YRI","chrAll.local.NAM","chrX.local.NAM")
boxplot(mxl[mxl$unrelated,cols],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"),cex.lab=0.85,
        main="Proportion Avg Local Ancestry in 53 Unrelated MXL Subjects\nFor Autosomes and X Chr Separately")
dev.off()


#### make "manhattan" plot of admixture res by chromosome
euro <- seq(from=79,to=145,by=3)
names(mxl[euro]) # good!

afr <- euro+1
nam <- afr+1

names(mxl[afr]); names(mxl[nam])
# so in order from chr 1-22, then x

# will have different order of samples depending on ancestry for each chr
# only plot the 53 unrelated

toPlot <- matrix(NA,nrow=3,ncol=53*23)
for(i in 1:23){
  ord <- order(mxl[mxl$unrelated,euro[i]])
  er <- mxl[mxl$unrelated,euro[i]][ord]
  af <- mxl[mxl$unrelated,euro[i]+1][ord]
  nm <- mxl[mxl$unrelated,euro[i]+2][ord]
  toPlot[,(i*53-52):(i*53)] <- rbind(er,af,nm)
}

pdf("../frappe_allChrs_manh_localAncest.pdf",width=35)
barplot(height=toPlot,col=c("blue","red","green"),space=0,axes=F,border=F)
axis(2)
axis(1,at=seq(from=27,to=ncol(toPlot),by=53),lab=paste("Chr",c(1:22,"X")))
dev.off()

## do just autosomes and x chr next to eachother
# take results for chrAll and Xchr
# cols 73-75, 70-72 are xchr
pdf("../frappe_auto_xChr_localAncest.pdf",width=14)
par(mfrow=c(1,2))
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]
# autosomes
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap MXL Avg Local Autosomal Ancestry")
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap MXL Avg Local X Chromsome Ancestry")
axis(2)
dev.off()

##
### make a plot of autos and x chr for mxl and asw in 2 panels
pdf("../frappe_auto_xChr_aswMXL_localAncest.pdf")
mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]
# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)
toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA)
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)
axis(1,at=c(26,79),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)


mxl <- get(load("asw_estimates_withLocal.RData"))
dim(mxl) # 87 148
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]

# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)

toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA)
axis(2)
axis(1,at=c(22,67),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("B", side=3, line=.75, adj=0, cex=1.3)
dev.off()


## add plot titles 
pdf("../frappe_auto_xChr_aswMXL_localAncest_withTitles.pdf")
mxl <- get(load("asw_estimates_withLocal.RData"))
par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]
# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)
toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA,main="HapMap ASW Estimated Ancestry")
axis(2)
axis(1,at=c(26,79),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)


mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
dim(mxl) # 87 150
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]

# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)

toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA,main="HapMap MXL Estimated Ancestry")
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)
axis(1,at=c(22,67),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("B", side=3, line=.75, adj=0, cex=1.3)
dev.off()


## make a png version of this

## add plot titles 
png("../frappe_auto_xChr_aswMXL_localAncest_withTitles.png")
mxl <- get(load("asw_estimates_withLocal.RData"))
par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]
# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)
toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA,main="HapMap ASW Estimated Ancestry")
axis(2)
axis(1,at=c(26,79),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)


mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
dim(mxl) # 87 150
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]

# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)

toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA,main="HapMap MXL Estimated Ancestry")
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)
axis(1,at=c(22,67),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("B", side=3, line=.75, adj=0, cex=1.3)
dev.off()

## make plot of only MXL samples 
pdf("../frappe_auto_xChr_MXL_localAncest_withTitles.pdf",width=11)
mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
dim(mxl) # 87 150
toPlOrd <- order(mxl[mxl$unrelated,"chrAll.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrAll.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrAll.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrAll.local.NAM"][toPlOrd]

# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,"chrX.local.CEU"])
toPl1 <- mxl[mxl$unrelated,"chrX.local.CEU"][toPlOrd]
toPl2 <- mxl[mxl$unrelated,"chrX.local.YRI"][toPlOrd]
toPl3 <- mxl[mxl$unrelated,"chrX.local.NAM"][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)

toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA,main="HapMap MXL Average Local Ancestry",
        cex.main=1.3)
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2,cex.axis=1.3)
axis(1,at=c(25,82),lab=c("Autosomal","X Chromosome"),tick=F,cex.axis=1.3)
dev.off()



mxl <- get(load("mxl_estimates_withLocal_v2.RData"))

range(mxl[mxl$unrelated,"chrX.local.CEU"]); sd(mxl[mxl$unrelated,"chrX.local.CEU"]) # euro
# 0 1
# 0.2407127
range(mxl[mxl$unrelated,"chrX.local.YRI"]); sd(mxl[mxl$unrelated,"chrX.local.YRI"]) # afr
# 0.000000 0.296482
# 0.05962235
range(mxl[mxl$unrelated,"chrX.local.NAM"]); sd(mxl[mxl$unrelated,"chrX.local.NAM"]) # nAm
# 0 1
# 0.2389882

mean(mxl[mxl$unrelated,"chrX.local.YRI"]) # 0.04520005
mean(mxl[mxl$unrelated,"chrX.local.NAM"]) # 0.5812864


####### try simple t-test
# pool all autosomal ancestries together, compare w all x chr ancestries
# a simple t-test on the ancestries together, NOT paired


# Native American
x <- mxl[mxl$unrelated,"chrAll.local.NAM"]
y <- mxl[mxl$unrelated,"chrX.local.NAM"]
t <- t.test(x,y)
t$p.value # 0.001046136

# African
x <- mxl[mxl$unrelated,"chrAll.local.YRI"]
y <- mxl[mxl$unrelated,"chrX.local.YRI"]
t <- t.test(x,y)
t$p.value # 0.7062355

# European
x <- mxl[mxl$unrelated,"chrAll.local.CEU"]
y <- mxl[mxl$unrelated,"chrX.local.CEU"]
t <- t.test(x,y)
t$p.value # 0.001457832


#####
# try a parallel coords plot

mxlDF <- data.frame(rbind(mxl[,c(1:3)],mxl[,c(1:3)],mxl[,c(1:3)]))
mxlDF$ancest <- c(rep("CEU",nrow(mxl)),rep("YRI",nrow(mxl)),rep("NAM",nrow(mxl)))
mxlDF$chr1 <- NA; mxlDF$chr2 <- NA; mxlDF$chr3 <- NA
mxlDF$chr4 <- NA; mxlDF$chr5 <- NA; mxlDF$chr6 <- NA
mxlDF$chr7 <- NA; mxlDF$chr8 <- NA; mxlDF$chr9 <- NA
mxlDF$chr10 <- NA; mxlDF$chr11 <- NA; mxlDF$chr12 <- NA
mxlDF$chr13 <- NA; mxlDF$chr14 <- NA; mxlDF$chr15 <- NA
mxlDF$chr16 <- NA; mxlDF$chr17 <- NA; mxlDF$chr18 <- NA
mxlDF$chr19 <- NA; mxlDF$chr20 <- NA; mxlDF$chr21 <- NA
mxlDF$chr22 <- NA; mxlDF$chr23 <- NA
ct <- 5
for(i in 1:22){
  cols <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),
            paste("chr",i,".local.NAM",sep=""))
  mxlDF[,ct] <- c(mxl[,cols[1]],mxl[,cols[2]],mxl[,cols[3]])
  ct <- ct+1
}
mxlDF$chr23 <- c(mxl[,"chrX.local.CEU"],mxl[,"chrX.local.YRI"],mxl[,"chrX.local.NAM"])

plot(1:23,1:23,ylim=c(0,1),xlab="Chromosome", ylab="Proportion Ancestry",type="n")

cols <- seq(from=79,to=145,by=1)
colnames(mxl)[cols]
mxlDF$color <- c(rep("cyan",nrow(mxl)),rep("magenta",nrow(mxl)),rep("purple",nrow(mxl)))
pdf("../parCoord_localAncest_mxl.pdf",width=14,height=10)
parcoord(mxlDF[,5:(ncol(mxlDF)-1)],col=mxlDF$color)
dev.off()


library(scatterplot3d)
mxl3d <- mxl[,c(1:3,cols)]
scatterplot3d(1:nrow(mxl),1:23,mxl[mxl$unrelated,cols], main="3D Scatterplot")

rm(list=ls())

#####
## get empirical distribution

setwd("~/Documents/tim_stuff/r_code")
source("emp_fisher.R")
library(CAnD)

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
mxl <- mxl[mxl$unrelated,]
dim(mxl) # 53 150

niters <- 1e6
ceu <- emp_fisher(mxl[,seq(from=79,to=146,by=3)],niters)
nam <- emp_fisher(mxl[,seq(from=81,to=147,by=3)],niters)
yri <- emp_fisher(mxl[,seq(from=80,to=147,by=3)],niters)

# get the actual estimates for the dat
true_ceu <- CAnD(mxl[,seq(from=79,to=146,by=3)],bonfCorr=FALSE)
true_nam <- CAnD(mxl[,seq(from=81,to=147,by=3)],bonfCorr=FALSE)
true_yri <- CAnD(mxl[,seq(from=80,to=147,by=3)],bonfCorr=FALSE)

summary(ceu)
summary(nam)
summary(yri)

write.table(ceu,file="ceu_1Miters.txt",quote=FALSE,row.names=FALSE)
write.table(nam,file="nam_1Miters.txt",quote=FALSE,row.names=FALSE)
write.table(yri,file="yri_1Miters.txt",quote=FALSE,row.names=FALSE)

rm(list=ls())

####
# plot empirical distribution
# get pvalue from empirical distribution

setwd("/home/staff/mchughc/admixed_haplotypes")
library(readr); library(dplyr); library(tidyr)
library(CAnD); library(ggplot2)

rowmns <- function(x) {
  rowMeans(chrAncest[-c(x)])
}

pairedTtest <- function(x) {
  t.test(chrAncest[, x], diff_means[, x], paired = TRUE)$p.value
}

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
mxl <- mxl[mxl$unrelated,]
dim(mxl) # 53 150

numChrs <- 23
n <- 53

chrAncest <- mxl[,seq(from=79,to=146,by=3)]
diff_means <- vapply(seq_len(numChrs), rowmns, rep(0, n))
pval <- vapply(seq_len(numChrs), pairedTtest, 0)
true_ceu <- -2 * sum(log(pval))

chrAncest <- mxl[,seq(from=81,to=147,by=3)]
diff_means <- vapply(seq_len(numChrs), rowmns, rep(0, n))
pval <- vapply(seq_len(numChrs), pairedTtest, 0)
true_nam <- -2 * sum(log(pval))

chrAncest <- mxl[,seq(from=80,to=147,by=3)]
diff_means <- vapply(seq_len(numChrs), rowmns, rep(0, n))
pval <- vapply(seq_len(numChrs), pairedTtest, 0)
true_yri <- -2 * sum(log(pval))

emp_yri <- read_delim("yri_1Miters.txt",delim=" ")
emp_ceu <- read_delim("ceu_1Miters.txt",delim=" ")
emp_nam <- read_delim("nam_1Miters.txt",delim=" ")

# calculate empirical pvalue
(calc_yri <- sum(emp_yri$x>=true_yri)/nrow(emp_yri)) # 0.463362
(calc_ceu <- sum(emp_ceu$x>=true_ceu)/nrow(emp_ceu)) # 0.441221
(calc_nam <- sum(emp_nam$x>=true_nam)/nrow(emp_nam)) # 0.438566


pdf("empiricalDist_zscore_yri_1Miters.pdf")
ggplot(emp_yri,aes(x=x)) + geom_density() + geom_vline(x=true_yri,color="red") + theme_bw() +
  xlab("Fishers Statistic for 1M Bootstrap Samples") + 
  ggtitle("YRI Null Empirical Distribution of Fisher's CAnD Statistic") +
  annotate("text",x=87,y=0.047,label=paste0("observed statistic = ",format(true_yri,digits=3))) + 
  annotate("text",x=86,y=0.045,label=paste0("empirical p-value = ",format(calc_yri,digits=4)))
dev.off()

pdf("empiricalDist_zscore_ceu_1Miters.pdf")
ggplot(emp_ceu,aes(x=x)) + geom_density() + geom_vline(x=true_ceu,color="red") + theme_bw() +
  xlab("Fishers Statistic for 1M Bootstrap Samples") + 
  ggtitle("CEU Null Empirical Distribution of Fisher's CAnD Statistic") +
  annotate("text",x=250,y=0.01,label=paste0("observed statistic = ",format(true_ceu,digits=3))) + 
  annotate("text",x=250,y=0.011,label=paste0("empirical p-value = ",format(calc_ceu,digits=4)))
dev.off()

pdf("empiricalDist_zscore_nam_1Miters.pdf")
ggplot(emp_nam,aes(x=x)) + geom_density() + geom_vline(x=true_nam,color="red") + theme_bw() +
  xlab("Fishers Statistic for 1M Bootstrap Samples") + 
  ggtitle("NAM Null Empirical Distribution of Fisher's CAnD Statistic") +
  annotate("text",x=250,y=0.01,label=paste0("observed statistic = ",format(true_nam,digits=3))) + 
  annotate("text",x=250,y=0.011,label=paste0("empirical p-value = ",format(calc_nam,digits=4)))
dev.off()

# wow, so empirical pvalue is nothing.
# do for like 1K iterations and save the column resamples too so we can see what's happening w chr 23

# do for NAM
ancestProp <- mxl[,seq(from=81,to=147,by=3)]
niters=1000
stats <- rep(NA,niters)
chrSamp <- matrix(NA,nrow=23,ncol=niters)
for(i in 1:niters){
  
  chr <- sample(1:23,23,replace=TRUE)
  #    chr <- chrIter[,i]
  chrAncest <- ancestProp[,chr]
  
  diff_means <- vapply(seq_len(numChrs), rowmns, rep(0, n))
  pval <- vapply(seq_len(numChrs), pairedTtest, 0)
  
  # want to store cand_stat = -2*sum(log(pval))
  stats[i] <- -2 * sum(log(pval))
  chrSamp[,i] <- chr
}

num23s <- apply(chrSamp,2,function(x){sum(x==23)})
table(num23s) 
#  0   1   2   3   4   5 
#359 387 192  51   9   2 

# look at the stats for each of these times
summary(stats[num23s==0])
summary(stats[num23s==1])
summary(stats[num23s==2])
summary(stats[num23s==3])
summary(stats[num23s==4])
summary(stats[num23s==5])
# yep, these correlate to the bumps in the distribution...

rm(list=ls())


#####
# try a new method of combining p-values
# want to estimate the correlation of our t-statistics

library(corpcor)

calc_combP <- function(dat){
  sig2.i <- apply(dat,1,var) # this is individual level variance
  
  # get the w_cc',i for each individual i
  m <- ncol(dat)
  tmp1 <- expand.grid(1:m,1:m) # all pairwise combos of chrs
  # remove rows where chrs are =
  cc <- tmp1[,1]==tmp1[,2]
  tmp1 <- tmp1[!cc,]
  w <- rep(NA,nrow(dat))
  for(i in 1:nrow(dat)){
    abar <- mean(as.numeric(dat[i,]))
    tmp2 <- cbind(as.numeric(dat[i,tmp1[,1]])-abar,as.numeric(dat[i,tmp1[,2]])-abar)
    w[i] <- mean(tmp2[,1]*tmp2[,2])
  }
  # w_{cc',i} should be zero under the null, since we are adj for the mean ancestry w/in an individ
  # this is relative to an individ, so there will be no correlation
  # if this is relative to the population itself, there will be correlation
  
  # need a sig2_c for each pair of chromosomes
  sig2.c <- rep(NA,ncol(dat))
  tstat <- rep(NA,ncol(dat))
  for(i in 1:m){
    tmp <- dat[,-i]
    pool <- apply(tmp,1,mean)
    d <- dat[,i]-pool
    dbar <- mean(d)
    #  sdd <- sd(d)/sqrt(n)
    #  tstat[i] <- dbar/sdd
    
    sig2.c[i] <- sum((d-dbar)^2)/(n*(n-1))
    tstat[i] <- dbar/sqrt(sig2.c[i])
    
    #  tstat[i] <- t.test(nam[,i],pool,paired=TRUE)$statistic
  }
  
  m <- ncol(dat)
  sig.matrix <- diag(nrow=m,ncol=m)
  for(i in 1:m){
    for(j in 1:m){
      denom <- n^2*sqrt(sig2.c[i])*sqrt(sig2.c[j])
      sig.matrix[i,j] <- (1/denom)*sum((m/(m-1)^2)*(w-sig2.i))
    }
  }
  diag(sig.matrix) <- 1
  mean(sig.matrix[lower.tri(sig.matrix)]) # -0.054286
  
  # calculate new stat
  (newstat <- tstat%*%solve(sig.matrix)%*%tstat) # 49.37731
  return(pchisq(newstat,df=m,lower.tail=FALSE))
}

# calculate the pairwise correlation between ancestry estimates for each individual
setwd("~/Documents/tim_stuff/r_code")
mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
mxl <- mxl[mxl$unrelated,]
dim(mxl) # 53 150

n <- nrow(mxl)
m <- ncol(mxl)

nam <- mxl[,seq(from=81,to=147,by=3)]
colnames(nam) # ok, good - chrs 1-22+X
dim(nam) # 53 23

(p <- calc_combP(nam)) # 0.001111169

namAuto <- nam[,1:22]
calc_combP(namAuto) # 0.4647446

###
eur <- mxl[,seq(from=79,to=145,by=3)]
colnames(eur) # ok, good - chrs 1-22+X
dim(eur) # 53 23

calc_combP(eur) # 0.00115851

eurAuto <- eur[,1:22]
calc_combP(eurAuto) # 0.5999746
# this is for EUR genome-wide


## afr now
afr <- mxl[,seq(from=80,to=147,by=3)]
colnames(afr) # ok, good - chrs 1-22+X
dim(afr)# 53 23

calc_combP(afr) # 0.9145613

afrAuto <- afr[,1:22]
calc_combP(afrAuto) # 0.9077587
# this is for AFR genome-wide


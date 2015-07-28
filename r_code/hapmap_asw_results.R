# hapmap ASW ancestry results

## load in all data
asw <- get(load("hapmap_aswOnly_estimates.RData"))
dim(asw) # 87 76


# check proportion native american --
png("../ancestry_differences/plots/asw_chrall_chrX_hgdp.png")
plot(asw[,"chrAll.hgdp"],asw[,"chrX.hgdp"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of HGDP Americas Ancestry")
points(asw[asw$unrelated,"chrAll.hgdp"],asw[asw$unrelated,"chrX.hgdp"],col="red")
points(asw[!asw$unrelated,"chrAll.hgdp"],asw[!asw$unrelated,"chrX.hgdp"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
dev.off()

# make a barplot of the x chr results
png("../ancestry_differences/plots/asw_chrx_frappe.png")
xchr <- asw[,c("chrX.ceu","chrX.yri","chrX.hgdp")]
xchr <- xchr[order(xchr$chrX.ceu,xchr$chrX.hgdp,xchr$chrX.yri),]
tt <- t(xchr)
colnames(tt) <- rep("",ncol(tt))
barplot(tt,col=c("blue","red","green"),xlab="Individual",main="HapMap ASW Estimated Ancestry\nX Chromosome")
dev.off()

l <- lm(asw[asw$unrelated,"chrX.hgdp"]~asw[asw$unrelated,"chrAll.hgdp"])
png("../ancestry_differences/plots/asw_chrall_chrX_hgdp_unrel.png")
plot(asw[,"chrAll.hgdp"],asw[,"chrX.hgdp"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of HGDP Americas Ancestry\nFor 45 Unrelated HapMap ASW")
points(asw[asw$unrelated,"chrAll.hgdp"],asw[asw$unrelated,"chrX.hgdp"],col="black")
#points(mxl[!mxl$unrelated,"chrAll.hgdp"],mxl[!mxl$unrelated,"chrX.hgdp"],col="black")
#legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
abline(summary(l)$coef[1,1],summary(l)$coef[2,1])
dev.off()

plot(abs(asw[,"chrAll.hgdp"]-asw[,"chrX.hgdp"]),
     type="n",main="Proportion of HGDP Americas Ancestry\nFor 45 Unrelated HapMap ASW")
points(abs(asw[asw$unrelated,"chrAll.hgdp"]-asw[asw$unrelated,"chrX.hgdp"]),col="black")

png("../ancestry_differences/plots/asw_chr15_chr8_hgdp.png")
plot(asw[,"chr15.hgdp"],asw[,"chr8.hgdp"],xlab="Chromosome 15",ylab="Chromosome 8",
     type="n",main="Proportion of HGDP Americas Ancestry")
points(asw[asw$unrelated,"chr15.hgdp"],asw[asw$unrelated,"chr8.hgdp"],col="red")
points(asw[!asw$unrelated,"chr15.hgdp"],asw[!asw$unrelated,"chr8.hgdp"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1)
dev.off()

png("../ancestry_differences/plots/asw_chr15_chrX_hgdp.png")
plot(asw[,"chr15.hgdp"],asw[,"chrX.hgdp"],xlab="Chromosome 15",ylab="X Chromosome",
     type="n",main="Proportion of HGDP Americas Ancestry")
points(asw[asw$unrelated,"chr15.hgdp"],asw[asw$unrelated,"chrX.hgdp"],col="red")
points(asw[!asw$unrelated,"chr15.hgdp"],asw[!asw$unrelated,"chrX.hgdp"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1)
dev.off()


var(asw$chrX.hgdp[asw$unrelated]) # 0.00128
var(asw$chrAll.hgdp[asw$unrelated]) # 0.001000
t.test(asw$chrX.hgdp[asw$unrelated],asw$chrAll.hgdp[asw$unrelated])
#t = -0.0931, df = 86.733, p-value = 0.9261
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.01479572  0.01347215
#sample estimates:
#  mean of x  mean of y 
#0.02187775 0.02253954 

var(asw$chr8.hgdp[asw$unrelated]) # 0.002596
var(asw$chr15.hgdp[asw$unrelated]) # 0.0031869
t.test(asw$chr8.hgdp[asw$unrelated],asw$chr15.hgdp[asw$unrelated])
#t = -0.3246, df = 87.091, p-value = 0.7463
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.02621073  0.01885227
#sample estimates:
#  mean of x  mean of y 
#0.02298864 0.02666786 


##################
# do analysis comparing each chr to the pool of all other chrs
res <- data.frame(matrix(NA,nrow=23,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=6,to=ncol(asw)-4,by=3)]
for(i in 1:23){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*23)
}       
res # all bonf p-vals are >1!

resEuro <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=4,to=ncol(asw)-4,by=3)]       
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*23)
  }       
resEuro # still nothing

resAfr <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=5,to=ncol(asw)-4,by=3)]       
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resAfr[i,] <- c(i,t$p.value,t$p.value*23)
}       
resAfr # nothing!

library(xtable)

resAfr$bonf_pvalue[resAfr$bonf_pvalue>1]<-1
resEuro$bonf_pvalue[resEuro$bonf_pvalue>1]<-1
res$bonf_pvalue[res$bonf_pvalue>1]<-1
xtable(cbind(res[,2:3],resEuro[,2:3],resAfr[,2:3]),digits=3)


#png("paired_ttest_pools.png")
pdf("../ancestry_differences/plots/asw_paired_ttest_pools.pdf")
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
pdf("../ancestry_differences/plots/asw_paired_ttest_bonfCorr_pools.pdf")
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
pdf("../ancestry_differences/plots/asw_paired_ttest_bonfAndNot.pdf")
par(mfrow=c(2,1))
plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="P-values from Paired T-Tests")
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)
plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="Bonferroni Corrected P-values from Paired T-Tests")
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
dev.off()



#### make boxplots of ancestry proportions for each ancestry subpop
pdf("../ancestry_differences/plots/asw_boxplot_ancestryProp.pdf",width=10)
boxplot(asw[asw$unrelated,c(73,70,74,71,75,72)],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"),cex.lab=0.85,
        main="Proportion Ancestry in 45 Unrelated ASW Subjects\nFor Autosomes and X Chr Separately")
dev.off()

## make boxplots of both ASW and MXL samples together
# make sure x-axis labels all show up
pdf("../ancestry_differences/plots/asw_mxl_boxplot_ancestryProp.pdf",width=11,height=10)
par( mai=c(0.5, 0.65, 0.4, 0.15), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))
boxplot(asw[asw$unrelated,c(73,70,74,71,75,72)],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"))
mtext("A", side=3, line=0.75,adj=0,cex=1.3)
mxl <- get(load("hapmap_mxlOnly_estimates.RData"))
dim(mxl) # 86 78
boxplot(mxl[mxl$unrelated,c(73,70,74,71,75,72)],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"))
mtext("B", side=3, line=0.75,adj=0,cex=1.3)
dev.off()



#### make "manhattan" plot of admixture res by chromosome
euro <- seq(from=4,to=71,by=3)
names(asw[euro]) # good!

afr <- euro+1
nam <- afr+1

names(asw[afr]); names(asw[nam])
# so in order from chr 1-22, then x

# will have different order of samples depending on ancestry for each chr
# only plot the 45 unrelated

## do just autosomes and x chr next to eachother
# take results for chrAll and Xchr
# cols 73-75, 70-72 are xchr
pdf("../asw_frappe_auto_xChr.pdf",width=14)
par(mfrow=c(1,2))
toPlOrd <- order(asw[asw$unrelated,73])
toPl1 <- asw[asw$unrelated,73][toPlOrd]
toPl2 <- asw[asw$unrelated,74][toPlOrd]
toPl3 <- asw[asw$unrelated,75][toPlOrd]
# autosomes
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap ASW Autosomal Ancestry")
legend("left",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)

#xchr
toPlOrd <- order(asw[asw$unrelated,70])
toPl1 <- asw[asw$unrelated,70][toPlOrd]
toPl2 <- asw[asw$unrelated,71][toPlOrd]
toPl3 <- asw[asw$unrelated,72][toPlOrd]
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap ASW X Chromsome Ancestry")
axis(2)
dev.off()

# x chr summary statistics
range(asw[asw$unrelated,70]); sd(asw[asw$unrelated,70]) # euro
#[1] 1.11722e-12 6.65364e-01
#[1] 0.1435466
range(asw[asw$unrelated,71]); sd(asw[asw$unrelated,71]) # afr
#[1] 0.257434 0.999999
#[1] 0.1528913
range(asw[asw$unrelated,72]); sd(asw[asw$unrelated,72]) # nAm
#[1] 2.32179e-26 1.89480e-01
#[1] 0.03570911

# autosomal summary statistics
range(asw[asw$unrelated,73]); sd(asw[asw$unrelated,73]) # euro
# 0.06313174 0.39098417
# [1] 0.07488461
range(asw[asw$unrelated,74]); sd(asw[asw$unrelated,74]) # afr
# 0.5849860 0.9193735
# [1] 0.08082176
range(asw[asw$unrelated,75]); sd(asw[asw$unrelated,75]) # nAm
# 0.006093686 0.216595391
# [1] 0.03162537



mean(asw[asw$unrelated,72]) # 0.0219
mean(asw[asw$unrelated,71]) # 0.8238


####### try simple t-test
# pool all autosomal ancestries together, compare w all x chr ancestries
# a simple t-test on the ancestries together, NOT paired


# Native American
x <- asw[asw$unrelated,"chrAll.hgdp"]
y <- asw[asw$unrelated,"chrX.hgdp"]
t <- t.test(x,y)
t$p.value # 0.9260633

# African
x <- asw[asw$unrelated,"chrAll.yri"]
y <- asw[asw$unrelated,"chrX.yri"]
t <- t.test(x,y)
t$p.value # 0.2376916

# European
x <- asw[asw$unrelated,"chrAll.ceu"]
y <- asw[asw$unrelated,"chrX.ceu"]
t <- t.test(x,y)
t$p.value # 0.2174325

rm(list=ls())


#####
# Remake ttest results without plot titles, in one panel

asw <- get(load("hapmap_aswOnly_estimates.RData"))
dim(asw) # 87 76

# do analysis comparing each chr to the pool of all other chrs
res <- data.frame(matrix(NA,nrow=23,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=6,to=ncol(asw)-4,by=3)]
for(i in 1:23){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*23)
}       
res # all bonf p-vals are >1!

resEuro <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=4,to=ncol(asw)-4,by=3)]       
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*23)
}       
resEuro # still nothing

resAfr <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=5,to=ncol(asw)-4,by=3)]       
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resAfr[i,] <- c(i,t$p.value,t$p.value*23)
}       
resAfr # nothing!

res$bonf_pvalue[res$bonf_pvalue>1] <- 1
resAfr$bonf_pvalue[resAfr$bonf_pvalue>1] <- 1
resEuro$bonf_pvalue[resEuro$bonf_pvalue>1] <- 1

pdf("../ancestry_differences/plots/asw_paired_ttest_bonfAndNot_bothPops.pdf")
par( mai=c(0.5, 0.65, 0.4, 0.15), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,2))

plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(0,1.5))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)


plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(0,1.5))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
mtext("B", side=3, line=0.75,adj=0,cex=1.3)


### now add in MXL results
mxl <- get(load("hapmap_mxlOnly_estimates.RData"))
dim(mxl) # 86 78

res <- data.frame(matrix(NA,nrow=23,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=6,to=ncol(mxl)-2,by=3)]
for(i in 1:23){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*23)
}       
res # xchr bonf pvalue: 7.400872e-05

resEuro <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=4,to=ncol(mxl)-2,by=3)]       
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*23)
}       
resEuro # xchr bonf pvalue: 7.248653e-05

resAfr <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=5,to=ncol(mxl)-2,by=3)]       
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resAfr[i,] <- c(i,t$p.value,t$p.value*23)
}       
resAfr # nothing!

res$bonf_pvalue[res$bonf_pvalue>1] <- 1
resAfr$bonf_pvalue[resAfr$bonf_pvalue>1] <- 1
resEuro$bonf_pvalue[resEuro$bonf_pvalue>1] <- 1

plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(0,7))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:7,labels=c(1:7),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)

mtext("C", side=3, line=0.75,adj=0,cex=1.3)

plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(-0.01,5.8))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:6,labels=c(1:6),cex.axis=0.8)
mtext("D", side=3, line=0.75,adj=0,cex=1.3)

dev.off()






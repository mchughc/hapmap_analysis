# hapmap MXL ancestry results

## load in all data
mxl <- get(load("hapmap_mxlOnly_estimates.RData"))
dim(mxl) # 86 78


# check proportion native american --
png("mxl_chrall_chrX_hgdp.png")
plot(mxl[,"chrAll.hgdp"],mxl[,"chrX.hgdp"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of HGDP Americas Ancestry")
points(mxl[mxl$unrelated,"chrAll.hgdp"],mxl[mxl$unrelated,"chrX.hgdp"],col="red")
points(mxl[!mxl$unrelated,"chrAll.hgdp"],mxl[!mxl$unrelated,"chrX.hgdp"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
dev.off()

# make a barplot of the x chr results
png("mxl_chrx_frappe.png")
xchr <- mxl[,c("chrX.ceu","chrX.yri","chrX.hgdp")]
xchr <- xchr[order(xchr$chrX.ceu,xchr$chrX.hgdp,xchr$chrX.yri),]
tt <- t(xchr)
colnames(tt) <- rep("",ncol(tt))
barplot(tt,col=c("blue","red","green"),xlab="Individual",main="HapMap MXL Estimated Ancestry\nX Chromosome")
dev.off()

l <- lm(mxl[mxl$unrelated,"chrX.hgdp"]~mxl[mxl$unrelated,"chrAll.hgdp"])
png("mxl_chrall_chrX_hgdp_unrel.png")
plot(mxl[,"chrAll.hgdp"],mxl[,"chrX.hgdp"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of HGDP Americas Ancestry\nFor 53 Unrelated HapMap MXL")
points(mxl[mxl$unrelated,"chrAll.hgdp"],mxl[mxl$unrelated,"chrX.hgdp"],col="black")
#points(mxl[!mxl$unrelated,"chrAll.hgdp"],mxl[!mxl$unrelated,"chrX.hgdp"],col="black")
#legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
abline(summary(l)$coef[1,1],summary(l)$coef[2,1])
dev.off()

plot(abs(mxl[,"chrAll.hgdp"]-mxl[,"chrX.hgdp"]),
     type="n",main="Proportion of HGDP Americas Ancestry\nFor 53 Unrelated HapMap MXL")
points(abs(mxl[mxl$unrelated,"chrAll.hgdp"]-mxl[mxl$unrelated,"chrX.hgdp"]),col="black")

png("mxl_chr15_chr8_hgdp.png")
plot(mxl[,"chr15.hgdp"],mxl[,"chr8.hgdp"],xlab="Chromosome 15",ylab="Chromosome 8",
     type="n",main="Proportion of HGDP Americas Ancestry")
points(mxl[mxl$unrelated,"chr15.hgdp"],mxl[mxl$unrelated,"chr8.hgdp"],col="red")
points(mxl[!mxl$unrelated,"chr15.hgdp"],mxl[!mxl$unrelated,"chr8.hgdp"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1)
dev.off()

png("mxl_chr15_chrX_hgdp.png")
plot(mxl[,"chr15.hgdp"],mxl[,"chrX.hgdp"],xlab="Chromosome 15",ylab="X Chromosome",
     type="n",main="Proportion of HGDP Americas Ancestry")
points(mxl[mxl$unrelated,"chr15.hgdp"],mxl[mxl$unrelated,"chrX.hgdp"],col="red")
points(mxl[!mxl$unrelated,"chr15.hgdp"],mxl[!mxl$unrelated,"chrX.hgdp"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1)
dev.off()


var(mxl$chrX.hgdp[mxl$unrelated]) # 0.0615
var(mxl$chrAll.hgdp[mxl$unrelated]) # 0.0249
t.test(mxl$chrX.hgdp[mxl$unrelated],mxl$chrAll.hgdp[mxl$unrelated])
#t = 3.3003, df = 88.145, p-value = 0.001395
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#  0.05301465 0.21348513 
#sample estimates:
#  mean of x mean of y 
#0.5882694 0.4550195 

var(mxl$chr8.hgdp[mxl$unrelated]) # 0.03385
var(mxl$chr15.hgdp[mxl$unrelated]) # 0.035578
t.test(mxl$chr8.hgdp[mxl$unrelated],mxl$chr15.hgdp[mxl$unrelated])
#t = 0.4407, df = 103.936, p-value = 0.6603
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#  -0.05582331  0.08772593 
#sample estimates:
#  mean of x mean of y 
#0.4708911 0.4549398 

var(mxl$chrX.hgdp[mxl$unrelated]) # 0.0615
var(mxl$chr15.hgdp[mxl$unrelated]) # 0.035578
t.test(mxl$chrX.hgdp[mxl$unrelated],mxl$chr15.hgdp[mxl$unrelated])
#t = 3.115, df = 97.071, p-value = 0.002419
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#  0.04837961 0.21827964 
#sample estimates:
#  mean of x mean of y 
#0.5882694 0.4549398 

##################
# do analysis comparing each chr to the pool of all other chrs
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

library(xtable)
xtable(cbind(res[,2:3],resEuro[,2:3],resAfr[,2:3]))

resAfr$bonf_pvalue[resAfr$bonf_pvalue>1]<-1
resEuro$bonf_pvalue[resEuro$bonf_pvalue>1]<-1
res$bonf_pvalue[res$bonf_pvalue>1]<-1
xtable(cbind(res[,2:3],resEuro[,2:3],resAfr[,2:3]))

## also adjust for false discovery rate
res$fdr_pvalue <- p.adjust(res$pvalue,method="fdr")
resAfr$fdr_pvalue <- p.adjust(resAfr$pvalue,method="fdr")
resEuro$fdr_pvalue <- p.adjust(resEuro$pvalue,method="fdr")


#png("paired_ttest_pools.png")
pdf("paired_ttest_pools.pdf")
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
pdf("paired_ttest_bonfCorr_pools.pdf")
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
pdf("../ancestry_differences/plots/paired_ttest_bonfAndNot_paper.pdf")
par( mai=c(0.65, 0.65, 0.4, 0.18), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))

plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(0,7))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:7,labels=c(1:7),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)

mtext("A", side=3, line=0.75,adj=0,cex=1.3)

plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(0,5.8))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:6,labels=c(1:6),cex.axis=0.8)
mtext("B", side=3, line=0.75,adj=0,cex=1.3)

dev.off()


#################
# do analysis comparing each chr to the pool of all other chrs
# excluding the X chromosome
res <- data.frame(matrix(NA,nrow=22,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=6,to=ncol(mxl)-7,by=3)]
for(i in 1:22){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*22)
}       

resEuro <- data.frame(matrix(NA,nrow=22,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=4,to=ncol(mxl)-7,by=3)]       
for(i in 1:22){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*22)
}       

resAfr <- data.frame(matrix(NA,nrow=22,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- mxl[mxl$unrelated,seq(from=5,to=ncol(mxl)-7,by=3)]       
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
pdf("paired_ttest_pools_exclXchr.pdf")
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
pdf("paired_ttest_bonfCorr_pools_exclXchr.pdf")
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
pdf("paired_ttest_bonfAndNot_exclXchr.pdf")
par(mfrow=c(2,1))
bigMax <- max(c(-log10(res$pvalue),-log10(resEuro$pvalue),-log10(resAfr$pvalue)))
plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="P-values from Paired T-Tests",ylim=c(0,bigMax))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:22,labels=c(1:22),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
bigMax <- max(c(-log10(res$bonf_pvalue),-log10(resEuro$bonf_pvalue),-log10(resAfr$bonf_pvalue)))
plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,
     main="Bonferroni Corrected P-values from Paired T-Tests",ylim=c(0,2))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:22,labels=c(1:22),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topright",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)
dev.off()


#### make boxplots of ancestry proportions for each ancestry subpop
pdf("../ancestry_differences/plots/boxplot_ancestryProp.pdf",width=10)
boxplot(mxl[mxl$unrelated,c(73,70,74,71,75,72)],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"),cex.lab=0.85,
        main="Proportion Ancestry in 53 Unrelated MXL Subjects\nFor Autosomes and X Chr Separately")
dev.off()


#### make "manhattan" plot of admixture res by chromosome
euro <- seq(from=4,to=71,by=3)
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
  
pdf("../frappe_allChrs_manh.pdf",width=35)
barplot(height=toPlot,col=c("blue","red","green"),space=0,axes=F,border=F)
axis(2)
axis(1,at=seq(from=27,to=ncol(toPlot),by=53),lab=paste("Chr",c(1:22,"X")))
dev.off()

## do just autosomes and x chr next to eachother
# take results for chrAll and Xchr
# cols 73-75, 70-72 are xchr
pdf("../frappe_auto_xChr.pdf",width=14)
par(mfrow=c(1,2))
toPlOrd <- order(mxl[mxl$unrelated,73])
toPl1 <- mxl[mxl$unrelated,73][toPlOrd]
toPl2 <- mxl[mxl$unrelated,74][toPlOrd]
toPl3 <- mxl[mxl$unrelated,75][toPlOrd]
# autosomes
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap MXL Autosomal Ancestry")
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,70])
toPl1 <- mxl[mxl$unrelated,70][toPlOrd]
toPl2 <- mxl[mxl$unrelated,71][toPlOrd]
toPl3 <- mxl[mxl$unrelated,72][toPlOrd]
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap MXL X Chromsome Ancestry")
axis(2)
dev.off()

##
### make a plot of autos and x chr for mxl and asw in 2 panels
pdf("../frappe_auto_xChr_aswMXL.pdf")
mxl <- get(load("hapmap_mxlOnly_estimates.RData"))
par( mai=c(0.5, 0.65, 0.4, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))
toPlOrd <- order(mxl[mxl$unrelated,73])
toPl1 <- mxl[mxl$unrelated,73][toPlOrd]
toPl2 <- mxl[mxl$unrelated,74][toPlOrd]
toPl3 <- mxl[mxl$unrelated,75][toPlOrd]
# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

toPlOrd <- order(mxl[mxl$unrelated,70])
toPl1 <- mxl[mxl$unrelated,70][toPlOrd]
toPl2 <- mxl[mxl$unrelated,71][toPlOrd]
toPl3 <- mxl[mxl$unrelated,72][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)
toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA)
legend("topleft",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)
axis(1,at=c(26,79),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)


mxl <- get(load("hapmap_aswOnly_estimates.RData"))
dim(mxl)
toPlOrd <- order(mxl[mxl$unrelated,73])
toPl1 <- mxl[mxl$unrelated,73][toPlOrd]
toPl2 <- mxl[mxl$unrelated,74][toPlOrd]
toPl3 <- mxl[mxl$unrelated,75][toPlOrd]
# autosomes
auts <- rbind(toPl1,toPl2,toPl3)

#xchr
toPlOrd <- order(mxl[mxl$unrelated,70])
toPl1 <- mxl[mxl$unrelated,70][toPlOrd]
toPl2 <- mxl[mxl$unrelated,71][toPlOrd]
toPl3 <- mxl[mxl$unrelated,72][toPlOrd]
xc <- rbind(toPl1,toPl2,toPl3)

toPl <- cbind(auts,xc)
barplot(height=toPl,col=c("blue","red","green"),space=0,axes=F,border=NA)
axis(2)
axis(1,at=c(22,67),lab=c("Autosomal","X Chromosome"),tick=F)
mtext("B", side=3, line=.75, adj=0, cex=1.3)
dev.off()






range(mxl[mxl$unrelated,70]); sd(mxl[mxl$unrelated,70]) # euro
#[1] 9.65339e-07 9.99999e-01
#[1] 0.2450824
range(mxl[mxl$unrelated,71]); sd(mxl[mxl$unrelated,71]) # afr
#[1] 2.45546e-33 2.96325e-01
#[1] 0.06070492
range(mxl[mxl$unrelated,72]); sd(mxl[mxl$unrelated,72]) # nAm
#[1] 9.87247e-07 9.99999e-01
#[1] 0.2480331

mean(mxl[mxl$unrelated,72]) # 0.588
mean(mxl[mxl$unrelated,71]) # 0.0434


####### try simple t-test
# pool all autosomal ancestries together, compare w all x chr ancestries
# a simple t-test on the ancestries together, NOT paired


# Native American
x <- mxl[mxl$unrelated,"chrAll.hgdp"]
y <- mxl[mxl$unrelated,"chrX.hgdp"]
t <- t.test(x,y)
t$p.value # 0.001395115

# try this by hand to be sure i'm getting the right thing
x1bar <- mean(x); x2bar <- mean(y)
n1 <- length(x); n2 <- length(y)
s1sq <- (1/(n1-1))*sum((x-x1bar)^2)
s2sq <- (1/(n2-1))*sum((y-x2bar)^2)
tstat <- (x1bar-x2bar)/(sqrt(s1sq/n1+s2sq/n2))
t$statistic; tstat # great!

# and the df calculation:
t$parameter # df=88.14495
(s1sq+s2sq)^2*(n1-1)/(s1sq^2+s2sq^2) # 88.14495; great!

# African
x <- mxl[mxl$unrelated,"chrAll.yri"]
y <- mxl[mxl$unrelated,"chrX.yri"]
t <- t.test(x,y)
t$p.value # 0.242359

# European
x <- mxl[mxl$unrelated,"chrAll.ceu"]
y <- mxl[mxl$unrelated,"chrX.ceu"]
t <- t.test(x,y)
t$p.value # 0.002733698




## do parametric test using ASW samples with local ancestry estimates

# hapmap ASW ancestry results

## load in all data
asw <- get(load("asw_estimates_withLocal.RData"))
dim(asw) # 87 148


# check proportion native american --
png("../ancestry_differences/plots/asw_chrall_chrX_hgdp_localAncest.png")
plot(asw[,"chrAll.local.NAM"],asw[,"chrX.local.NAM"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of Avg Local HGDP Americas Ancestry")
points(asw[asw$unrelated,"chrAll.local.NAM"],asw[asw$unrelated,"chrX.local.NAM"],col="red")
points(asw[!asw$unrelated,"chrAll.local.NAM"],asw[!asw$unrelated,"chrX.local.NAM"],col="black")
legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
dev.off()

# make a barplot of the x chr results
png("../ancestry_differences/plots/asw_chrx_frappe_localAncest.png")
xchr <- asw[,c("chrX.local.CEU","chrX.local.YRI","chrX.local.NAM")]
xchr <- xchr[order(xchr$chrX.local.CEU,xchr$chrX.local.NAM,xchr$chrX.local.YRI),]
tt <- t(xchr)
colnames(tt) <- rep("",ncol(tt))
barplot(tt,col=c("blue","red","green"),xlab="Individual",main="HapMap ASW Estimated Avg Local Ancestry\nX Chromosome")
dev.off()

l <- lm(asw[asw$unrelated,"chrX.local.NAM"]~asw[asw$unrelated,"chrAll.local.NAM"])
png("../ancestry_differences/plots/asw_chrall_chrX_hgdp_unrel_localAncest.png")
plot(asw[,"chrAll.local.NAM"],asw[,"chrX.local.NAM"],xlab="All Autosomal SNPs",ylab="X Chromosome",
     type="n",main="Proportion of Avg Local HGDP Americas Ancestry\nFor 45 Unrelated HapMap ASW")
points(asw[asw$unrelated,"chrAll.local.NAM"],asw[asw$unrelated,"chrX.local.NAM"],col="black")
#points(mxl[!mxl$unrelated,"chrAll.local.NAM"],mxl[!mxl$unrelated,"chrX.local.NAM"],col="black")
#legend("bottomright",c("Unrelated"),col="red",pch=1)
abline(0,1,lty=2,col="grey")
abline(summary(l)$coef[1,1],summary(l)$coef[2,1])
dev.off()



var(asw$chrX.local.NAM[asw$unrelated]) # 0.001340183
var(asw$chrAll.local.NAM[asw$unrelated]) # 0.001024059
t.test(asw$chrX.local.NAM[asw$unrelated],asw$chrAll.local.NAM[asw$unrelated])
# t = 0.036, df = 86.454, p-value = 0.9714
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.01414720  0.01466916
# sample estimates:
#   mean of x  mean of y 
# 0.01470146 0.01444048 

var(asw$chr8.local.NAM[asw$unrelated]) # 0.002502811
var(asw$chr15.local.NAM[asw$unrelated]) # 0.002773682
t.test(asw$chr8.local.NAM[asw$unrelated],asw$chr15.local.NAM[asw$unrelated])
# t = -0.0546, df = 87.769, p-value = 0.9566
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02211096  0.02092917
# sample estimates:
#   mean of x  mean of y 
# 0.01558334 0.01617423 


##################
# do analysis comparing each chr to the pool of all other chrs
res <- data.frame(matrix(NA,nrow=23,ncol=3))
names(res) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=79,to=145,by=3)]
colnames(sm) # good -- all autosomal and x chr NAM estimates
for(i in 1:23){
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  res[i,] <- c(i,t$p.value,t$p.value*23)
}       
res # all bonf p-vals are >1!

resEuro <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resEuro) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=77,to=143,by=3)]       
colnames(sm)
for(i in 1:23){     
  tmp <- sm[,-i]
  pool <- apply(tmp,1,mean)
  t <- t.test(sm[,i],pool,paired=TRUE)
  resEuro[i,] <- c(i,t$p.value,t$p.value*23)
}       
resEuro # still nothing

resAfr <- data.frame(matrix(NA,nrow=23,ncol=3))
names(resAfr) <- c("chr","pvalue","bonf_pvalue")
sm <- asw[asw$unrelated,seq(from=78,to=144,by=3)]       
colnames(sm)
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
pdf("../ancestry_differences/plots/asw_paired_ttest_pools_localAncest.pdf")
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
pdf("../ancestry_differences/plots/asw_paired_ttest_bonfCorr_pools_localAncest.pdf")
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
pdf("../ancestry_differences/plots/asw_paired_ttest_bonfAndNot_localAncest.pdf")
par( mai=c(0.65, 0.65, 0.4, 0.15), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))

miny <- min(c(-log10(res$pvalue),-log10(resEuro$pvalue),-log10(resAfr$pvalue)))
maxy <- max(c(-log10(res$pvalue),-log10(resEuro$pvalue),-log10(resAfr$pvalue)))
plot(res$chr,-log10(res$pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(miny,2.2))
points(resEuro$chr,-log10(resEuro$pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
legend("topleft",c("European","African","Native American"),pch=c(2,3,19),
       col=c("blue","red","green"),cex=0.8)
mtext("A", side=3, line=0.75,adj=0,cex=1.3)

miny <- min(c(-log10(res$bonf_pvalue),-log10(resEuro$bonf_pvalue),-log10(resAfr$bonf_pvalue)))
maxy <- max(c(-log10(res$bonf_pvalue),-log10(resEuro$bonf_pvalue),-log10(resAfr$bonf_pvalue)))
plot(res$chr,-log10(res$bonf_pvalue),pch=19,col="green",xlab="Chromosome",axes=FALSE,
     ylab=expression(paste(-log[10],"(p-value)",sep="")),frame.plot=TRUE,ylim=c(miny,1.5))
points(resEuro$chr,-log10(resEuro$bonf_pvalue),pch=2,col="blue")
points(resAfr$chr,-log10(resAfr$bonf_pvalue),pch=3,col="red")
abline(h=-log10(0.05),col="gray",lty=2)
axis(1,at=1:23,labels=c(1:22,"X"),cex.axis=0.75)
axis(2,at=1:5,labels=c(1:5),cex.axis=0.8)
mtext("B", side=3, line=0.75,adj=0,cex=1.3)

dev.off()



#### make boxplots of ancestry proportions for each ancestry subpop
pdf("../ancestry_differences/plots/asw_boxplot_ancestryProp_localAncest.pdf",width=10)
cols <- c("chrAll.local.CEU","chrX.local.CEU","chrAll.local.YRI","chrX.local.YRI","chrAll.local.NAM","chrX.local.NAM")
boxplot(asw[asw$unrelated,cols],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"),cex.lab=0.85,
        main="Proportion Avg Local Ancestry in 45 Unrelated ASW Subjects\nFor Autosomes and X Chr Separately")
dev.off()

## make boxplots of both ASW and MXL samples together
# make sure x-axis labels all show up
pdf("../ancestry_differences/plots/asw_mxl_boxplot_ancestryProp_localAncest.pdf",width=11,height=10)
par( mai=c(0.5, 0.65, 0.4, 0.15), mgp=c(2, 0.5, 0), tck=-0.03 ,mfrow=c(2,1))
#par(mfrow=c(2,1))
cols <- c("chrAll.local.CEU","chrX.local.CEU","chrAll.local.YRI","chrX.local.YRI","chrAll.local.NAM","chrX.local.NAM")
boxplot(asw[asw$unrelated,cols],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"))
mtext("A", side=3, line=0.75,adj=0,cex=1.3)
mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
dim(mxl) # 86 78
boxplot(mxl[mxl$unrelated,cols],names=c("Euro, Auto","Euro, X Chr","Afr, Auto","Afr, X Chr",
                                                        "Native Am, Auto","Native Am, X Chr"),
        ylab="Proportion Ancestry",col=c("purple",NA,"cyan",NA,"orange",NA),
        border=c("black","purple","black","cyan","black","orange"))
mtext("B", side=3, line=0.75,adj=0,cex=1.3)
dev.off()



#### make "manhattan" plot of admixture res by chromosome
euro <- seq(from=77,to=143,by=3)
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
pdf("../asw_frappe_auto_xChr_localAncest.pdf",width=14)
par(mfrow=c(1,2))
toPlOrd <- order(asw[asw$unrelated,73])
toPl1 <- asw[asw$unrelated,73][toPlOrd]
toPl2 <- asw[asw$unrelated,74][toPlOrd]
toPl3 <- asw[asw$unrelated,75][toPlOrd]
# autosomes
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap ASW Avg Local Autosomal Ancestry")
legend("left",c("European","African","Native American"),lty=1,col=c("blue","red","green"),lwd=6,bg="white")
axis(2)

#xchr
toPlOrd <- order(asw[asw$unrelated,70])
toPl1 <- asw[asw$unrelated,70][toPlOrd]
toPl2 <- asw[asw$unrelated,71][toPlOrd]
toPl3 <- asw[asw$unrelated,72][toPlOrd]
barplot(height=rbind(toPl1,toPl2,toPl3),col=c("blue","red","green"),space=0,axes=F,border=T,xlab="Individual",
        main="HapMap ASW Avg Local X Chromsome Ancestry")
axis(2)
dev.off()

# x chr summary statistics
range(asw[asw$unrelated,"chrX.local.CEU"]); sd(asw[asw$unrelated,"chrX.local.CEU"]) # euro
# 0.0000000 0.6740545
# 0.1384112
range(asw[asw$unrelated,"chrX.local.YRI"]); sd(asw[asw$unrelated,"chrX.local.YRI"]) # afr
# 0.3259455 1.0000000
# 0.1397125
range(asw[asw$unrelated,"chrX.local.NAM"]); sd(asw[asw$unrelated,"chrX.local.NAM"]) # nAm
# 0.0000000 0.1985928
# 0.03660851

# autosomal summary statistics
range(asw[asw$unrelated,"chrAll.local.CEU"]); sd(asw[asw$unrelated,"chrAll.local.CEU"]) # euro
# 0.07085973 0.40145501
# 0.0774005
range(asw[asw$unrelated,"chrAll.local.YRI"]); sd(asw[asw$unrelated,"chrAll.local.YRI"]) # afr
# 0.5852573 0.9202229
# 0.0824855
range(asw[asw$unrelated,"chrAll.local.NAM"]); sd(asw[asw$unrelated,"chrAll.local.NAM"]) # nAm
# 0.0003276719 0.2089644814
# 0.03200092


mean(asw[asw$unrelated,"chrX.local.NAM"]) # 0.0147
mean(asw[asw$unrelated,"chrX.local.YRI"]) # 0.8316


####### try simple t-test
# pool all autosomal ancestries together, compare w all x chr ancestries
# a simple t-test on the ancestries together, NOT paired


# Native American
x <- asw[asw$unrelated,"chrAll.local.NAM"]
y <- asw[asw$unrelated,"chrX.local.NAM"]
t <- t.test(x,y)
t$p.value # 0.971361

# African
x <- asw[asw$unrelated,"chrAll.local.YRI"]
y <- asw[asw$unrelated,"chrX.local.YRI"]
t <- t.test(x,y)
t$p.value # 0.1748055

# European
x <- asw[asw$unrelated,"chrAll.local.CEU"]
y <- asw[asw$unrelated,"chrX.local.CEU"]
t <- t.test(x,y)
t$p.value # 0.16207


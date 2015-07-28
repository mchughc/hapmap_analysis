
### compare the ADMIXTURE results for each individual, by chromosome, to the average of the local ancestry
## only have for the hapmap MXL

## Contents:
# 1. Process files
# 2. Plot diff in ancestry by chr
# 3. Plot correlations by chr, by individ 
# 4. Compare avg local ancestry across genome with ADMIX.all estimates
# 5. Redo assortative mating analysis with local admixture estimates
# 6. Plot genome wide local and ADMIX diff for all individs
# 7. Look at var of local ancestry est across individs
# 8. Plot genome wide local and ADMIX diff for all individs on X chr
# 9. Get mean (SD) for estimates by chromosome



#####
# 1. Process files

mxl <- get(load("hapmap_mxlOnly_estimates.RData"))

# merge in sex info
trios <- read.table("../mex_trios.txt",header=T,as.is=T)
dim(trios); head(trios) # 72 7
parents <- trios[is.element(trios$father,0)&is.element(trios$mother,0),]
dim(parents) # 48 7
table(table(parents$family)) # 24 of 2

mates <- mxl[is.element(mxl$V2,parents$coriell.id),]
dim(mates) # 48 78

table(table(mates$superfamily))
#2  4  6 
#19  1  1

t <- table(mates$superfamily)
mates[is.element(mates$superfamily,names(t)[t==4]),]
# take NA19657, NA19658
mates <- mates[!is.element(mates[,"V2"],c("NA19785","NA19786")),]
dim(mates)
mates[is.element(mates$superfamily,names(t)[t==6]),]
# take the two that have superfamilyUnrel==TRUE
mates <- mates[!is.element(mates[,"V2"],c("NA19684","NA19661","NA19660","NA19685")),]
dim(mates) # 42 78; great!

table(table(mates[,1])) # 21 pairs
table(table(mates$superfamily)) # 21 pairs

# get average local ancestry for MXL individuals by chromosome
for(i in 1:22){
  ancest <- read.table(paste("../assortative_mating/HAPMAP_MEX_CHR_ANCESTRY/HAPMAP_MEX_CHR",i,".txt",sep=""),header=T)
  stopifnot(colnames(ancest)==c("FID","IID","CEU","YRI","NAM"))
  colnames(ancest) <- c("FID","IID",paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),
                        paste("chr",i,".local.NAM",sep=""))
  mxl <- merge(mxl,ancest[,2:5],by.x="V2",by.y="IID",all.x=TRUE)
}

# get x chr too
ancest <- read.table("../assortative_mating/HAPMAP_MEX_CHR_ANCESTRY/HAPMAP_MEX_CHR23.txt",header=T)
stopifnot(colnames(ancest)==c("FID","IID","CEU","YRI","NAM"))
colnames(ancest) <- c("FID","IID","chrX.local.CEU","chrX.local.YRI",
                      "chrX.local.NAM")
mxl <- merge(mxl,ancest[,2:5],by.x="V2",by.y="IID",all.x=TRUE)

save(mxl,file="mxl_estimates_withLocal.RData")


#####
# 2. Plot diff in ancestry by chr

# summarize the difference between the local and admixture estimates
# make a plot for each chromosome; x axis is individ and y axis is diff between estimates
# 3 lines: one per ancestry

x <- 1:nrow(mxl)

pdf("../chrs1_3_difference.pdf",height=18)
par(mfrow=c(3,1))
diff <- c(mxl$chr1.local.CEU-mxl$chr1.ceu,mxl$chr1.local.YRI-mxl$chr1.yri,mxl$chr1.local.NAM-mxl$chr1.hgdp)
plot(x,y=mxl$chr1.local.CEU-mxl$chr1.ceu,xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 1",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl$chr1.local.CEU-mxl$chr1.ceu,col="cyan",lwd=2,type="l")
points(x,mxl$chr1.local.YRI-mxl$chr1.yri,col="magenta",lwd=2,lty=2,type="l")
points(x,mxl$chr1.local.NAM-mxl$chr1.hgdp,col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topleft",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

diff <- c(mxl$chr2.local.CEU-mxl$chr2.ceu,mxl$chr2.local.YRI-mxl$chr2.yri,mxl$chr2.local.NAM-mxl$chr2.hgdp)
plot(x,y=mxl$chr2.local.CEU-mxl$chr2.ceu,xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 2",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl$chr2.local.CEU-mxl$chr2.ceu,col="cyan",lwd=2,type="l")
points(x,mxl$chr2.local.YRI-mxl$chr2.yri,col="magenta",lwd=2,lty=2,type="l")
points(x,mxl$chr2.local.NAM-mxl$chr2.hgdp,col="blue",lwd=2,lty=3,type="l")
abline(h=0)

diff <- c(mxl$chr3.local.CEU-mxl$chr3.ceu,mxl$chr3.local.YRI-mxl$chr3.yri,mxl$chr3.local.NAM-mxl$chr3.hgdp)
plot(x,y=mxl$chr3.local.CEU-mxl$chr3.ceu,xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 3",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl$chr3.local.CEU-mxl$chr3.ceu,col="cyan",lwd=2,type="l")
points(x,mxl$chr3.local.YRI-mxl$chr3.yri,col="magenta",lwd=2,lty=2,type="l")
points(x,mxl$chr3.local.NAM-mxl$chr3.hgdp,col="blue",lwd=2,lty=3,type="l")
abline(h=0)
dev.off()


pdf("../chrs4_6_difference.pdf",height=18)
par(mfrow=c(3,1))

i=4
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 4",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topleft",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

i=5
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 5",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=6
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 6",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

dev.off()


pdf("../chrs7_9_difference.pdf",height=18)
par(mfrow=c(3,1))

i=7
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 7",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topright",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

i=8
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 8",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=9
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 9",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

dev.off()


pdf("../chrs10_12_difference.pdf",height=18)
par(mfrow=c(3,1))

i=10
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 10",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topright",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

i=11
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 11",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=12
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 12",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

dev.off()


pdf("../chrs13_15_difference.pdf",height=18)
par(mfrow=c(3,1))

i=13
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 13",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topright",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

i=14
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 14",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=15
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 15",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

dev.off()



pdf("../chrs16_18_difference.pdf",height=18)
par(mfrow=c(3,1))

i=16
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 16",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topright",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

i=17
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 17",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=18
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 18",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

dev.off()


pdf("../chrs19_22_difference.pdf",height=18)
par(mfrow=c(4,1))

i=19
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 19",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topright",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)

i=20
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 20",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=21
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 21",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)

i=22
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="Chr 22",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
dev.off()


pdf("../chrX_difference.pdf")
i="X"
cols1 <- c(paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),paste("chr",i,".local.NAM",sep=""))
cols2 <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))

diff <- mxl[,cols1]-mxl[,cols2]
plot(x,y=mxl[,cols1[1]]-mxl[,cols2[1]],xlab="Individual",ylab="Mean Local Ancestry - ADMIXTURE Estimate",
     main="X Chromosome",type="n",ylim=c(min(diff),max(diff)))
points(x,mxl[,cols1[1]]-mxl[,cols2[1]],col="cyan",lwd=2,type="l")
points(x,mxl[,cols1[2]]-mxl[,cols2[2]],col="magenta",lwd=2,lty=2,type="l")
points(x,mxl[,cols1[3]]-mxl[,cols2[3]],col="blue",lwd=2,lty=3,type="l")
abline(h=0)
legend("topleft",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)
dev.off()


#### ok, so got the differences
# biggest differences are on chr 22 in a few individuals
# generally small differences except for one or two individuals per chromosome
# on the whole, admixture is under estimating the amount CEU, over estimating NAm and about equal on YRI


#####
# 3. Plot correlations by chr, by individ

corrAllIndivids <- data.frame(matrix(NA,nrow=22,ncol=3))
colnames(corrAllIndivids) <- c("CEU","YRI","NAm")

# get correlation values of the estimates, for each ancestry, by chromosome across all individs
for(i in 1:22){
  col.local <- paste("chr",i,".local.CEU",sep="")
  col.admix <- paste("chr",i,".ceu",sep="")
  corrAllIndivids[i,"CEU"] <- cor(mxl[,col.local],mxl[,col.admix])

  col.local <- paste("chr",i,".local.YRI",sep="")
  col.admix <- paste("chr",i,".yri",sep="")
  corrAllIndivids[i,"YRI"] <- cor(mxl[,col.local],mxl[,col.admix])
  
  col.local <- paste("chr",i,".local.NAM",sep="")
  col.admix <- paste("chr",i,".hgdp",sep="")
  corrAllIndivids[i,"NAm"] <- cor(mxl[,col.local],mxl[,col.admix])
}

corrAllChrs <- data.frame(matrix(NA,nrow=nrow(mxl),ncol=3))
colnames(corrAllChrs) <- c("CEU","YRI","NAm")

# get correlation values of the estimates, for each ancestry, by individ across the genome
for(i in 1:nrow(mxl)){
  col.local <- seq(from=79,to=144,by=3)
  col.admix <- seq(from=4,to=69,by=3)
  corrAllChrs[i,"CEU"] <- cor(as.numeric(mxl[i,col.local]),as.numeric(mxl[i,col.admix]))
  
  col.local <- col.local+1
  col.admix <- col.admix+1
  corrAllChrs[i,"YRI"] <- cor(as.numeric(mxl[i,col.local]),as.numeric(mxl[i,col.admix]))
  
  col.local <- col.local+1
  col.admix <- col.admix+1
  corrAllChrs[i,"NAm"] <- cor(as.numeric(mxl[i,col.local]),as.numeric(mxl[i,col.admix]))
}

x <- 1:nrow(corrAllIndivids)

pdf("../correlation_localvsAdmix_acrossIndivids.pdf")
plot(x,corrAllIndivids[,"CEU"],type="b",col="cyan",ylim=c(min(corrAllIndivids),max(corrAllIndivids)),
     lwd=2,xlab="Chromosome",ylab="Correlation of Ancestry Estimates Across All MXL Individs",pch=19)
points(x,corrAllIndivids[,"YRI"],type="b",col="magenta",lty=2,lwd=2,pch=19)
points(x,corrAllIndivids[,"NAm"],type="b",col="blue",lty=3,lwd=2,pch=19)
legend("bottomleft",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)
dev.off()

pdf("../correlation_localvsAdmix_acrossChrs.pdf")
x <- 1:nrow(corrAllChrs)
plot(x,corrAllChrs[,"CEU"],type="b",col="cyan",ylim=c(min(corrAllChrs),max(corrAllChrs)),
     lwd=2,xlab="Individual",ylab="Correlation of Ancestry Estimates Across All Chrs",pch=19)
points(x,corrAllChrs[,"YRI"],type="b",col="magenta",lty=2,lwd=2,pch=19)
points(x,corrAllChrs[,"NAm"],type="b",col="blue",lty=3,lwd=2,pch=19)
legend("bottomleft",c("CEU","YRI","NAm"),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)
dev.off()

rm(list=ls())


#####
# 4. Compare avg local ancestry across genome with ADMIX.all estimates

mxl <- get(load("mxl_estimates_withLocal.RData"))
col.local <- seq(from=79,to=144,by=3)
mxl$chrAll.local.CEU <- rowMeans(mxl[,col.local])

col.local <- col.local+1
mxl$chrAll.local.YRI <- rowMeans(mxl[,col.local])

col.local <- col.local+1
mxl$chrAll.local.NAM <- rowMeans(mxl[,col.local])

summary(mxl$chrAll.local.CEU-mxl$chrAll.ceu)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.045710 -0.001157  0.008376  0.010080  0.022100  0.059810 

summary(mxl$chrAll.local.YRI-mxl$chrAll.yri)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.026440 -0.008435 -0.003229 -0.003219  0.002353  0.023900 

summary(mxl$chrAll.local.NAM-mxl$chrAll.hgdp)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.067670 -0.017440 -0.005401 -0.006860  0.002681  0.039580 

save(mxl,file="mxl_estimates_withLocal_v2.RData")

rm(list=ls())


#####
# 5. Redo assortative mating analysis with local admixture estimates

source("getEst.R")

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
table(table(mxl[,2]))
#       1  2  3 
#       4  5 24
table(table(mxl$superfamily))
#1  2  3  6 10 
#3  5 19  1  1

# merge in sex info
trios <- read.table("../mex_trios.txt",header=T,as.is=T)
dim(trios); head(trios) # 72 7
parents <- trios[is.element(trios$father,0)&is.element(trios$mother,0),]
dim(parents) # 48 7
table(table(parents$family)) # 24 of 2

mates <- mxl[is.element(mxl$V2,parents$coriell.id),]
dim(mates) # 48 144

table(table(mates$superfamily))
#2  4  6 
#19  1  1

t <- table(mates$superfamily)
mates[is.element(mates$superfamily,names(t)[t==4]),]
# take NA19657, NA19658
mates <- mates[!is.element(mates[,1],c("NA19785","NA19786")),]
dim(mates) # 46 147
mates[is.element(mates$superfamily,names(t)[t==6]),]
# take the two that have superfamilyUnrel==TRUE
mates <- mates[!is.element(mates[,1],c("NA19684","NA19661","NA19660","NA19685")),]
dim(mates) # 42 147; great!

table(table(mates[,2])) # 21 pairs
table(table(mates$superfamily)) # 21 pairs

mates <- merge(mates,parents,by.y="coriell.id",by.x="V2")
dim(mates); head(mates) # 42 150
# now these are all the mate pairs, w sex info
sum(duplicated(mates$coriell.id)) # 0, so all unique subjects

table(mates$family,mates$sex) # good, so one F & one M per family
mates <- mates[order(mates$family),]

##
## do for HGDP/Native Am local ancestry estimates
nms <- paste("chr",c(1:22,"X","All"),".local.NAM",sep="")
ethnData <- mates[,c(nms,"V2","sex")]
names(ethnData) <- c(1:22,"x","all","coriell.id","sex")

res <- getEst(ethnData,5432,withX=TRUE,allChrs=TRUE)

chrxF <- ethnData[ethnData$sex=="F",c(1:22,"x","all")]
chrxM <- ethnData[ethnData$sex=="M",c(1:22,"x","all")]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:24){
  if(i==23){filen <- paste("../assortative_mating/plots/localAncestry_Corr_allchr_mateAncest.pdf",sep="")
            mainTitle <- paste("Correlation of Native Am local ancestry between mates\nacross all chromosomes")
            xlabAxis <- paste("Corr of All Chrs Native Am local ancestry")
            cx <- corrChr["all"]
  }else if(i==24){
    filen <- paste("../assortative_mating/plots/localAncestry_Corr_xchr_mateAncest.pdf",sep="")
    mainTitle <- paste("Correlation of Native Am local ancestry between mates\non X chromosomes")
    xlabAxis <- paste("Corr of X Chr Native Am local ancestry")
    cx <- corrChr["x"]
  }else{
    filen <- paste("../assortative_mating/plots/localAncestry_Corr_",i,"chr_mateAncest.pdf",sep="")
    mainTitle <- paste("Correlation of Native Am local ancestry between mates\non chromosome",i)
    xlabAxis <- paste("Corr of Chr",i,"Native Am local ancestry")
    cx <- corrChr[i]
  }
  pdf(filen)
  hist(res[i,],main=mainTitle,xlab=xlabAxis,
       sub="Empirical distribution of 5000 resamples",breaks=50)
  abline(v=cx,col="red",lwd=2)
  cx <- format(cx,digits=3)
  legend("topright",(paste("Obs correlation =",cx)),col="red",
         lty=1,bg="white")
  dev.off()
}

mean(res["x",]); sd(res["x",])# x chr results
# -0.0006108888 | 0.2217673
mean(res["all",]); sd(res["all",]) # all chrs results
#  0.002531879 | 0.2257373

(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.5003299
sum(res["x",]>=cx)/ncol(res) # 0.0132; there are 66 of 5000 greater
# test no random mating, ie assortative or disassortative
sum(abs(res["x",])>=cx)/ncol(res) # 0.0234; 117 of 5000 greater

(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.4384212
sum(res["all",]>=ca)/ncol(res) # 0.0224; there are 112 of 5000 greater
# test no random mating, ie assortative or disassortative
sum(abs(res["all",])>=ca)/ncol(res) # 0.0462; 231 of 5000 greater


###
## do for YRI
nms <- paste("chr",c(1:22,"X","All"),".local.YRI",sep="")
ethnData <- mates[,c(nms,"V2","sex")]
names(ethnData) <- c(1:22,"x","all","coriell.id","sex")

resAfr <- getEst(ethnData,5432,withX=TRUE,allChrs=TRUE)

chrxF <- ethnData[ethnData$sex=="F",c(1:22,"x","all")]
chrxM <- ethnData[ethnData$sex=="M",c(1:22,"x","all")]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:24){
  if(i==23){filen <- paste("../assortative_mating/plots/localAncestry_Corr_allchr_mateAncest_YRI.pdf",sep="")
            mainTitle <- paste("Correlation of African local ancestry between mates\nacross all chromosomes")
            xlabAxis <- paste("Corr of All Chrs African local ancestry")
            cx <- corrChr["all"]
  }else if(i==24){filen <- paste("../assortative_mating/plots/localAncestry_Corr_xchr_mateAncest_YRI.pdf",sep="")
                  mainTitle <- paste("Correlation of African local ancestry between mates\non X chromosome")
                  xlabAxis <- paste("Corr of X Chr African local ancestry")    
                  cx <- corrChr["x"]
  }else{
    filen <- paste("../assortative_mating/plots/localAncestry_Corr_",i,"chr_mateAncest_YRI.pdf",sep="")
    mainTitle <- paste("Correlation of African local ancestry between mates\non chromosome",i)
    xlabAxis <- paste("Corr of Chr",i,"African local ancestry")
    cx <- corrChr[i]
  }
  pdf(filen)
  hist(resAfr[i,],main=mainTitle,xlab=xlabAxis,
       sub="Empirical distribution of 5000 resamples",breaks=50)
  abline(v=cx,col="red",lwd=2)
  cx <- format(cx,digits=3)
  legend("topright",(paste("Obs correlation =",cx)),col="red",
         lty=1,bg="white")
  dev.off()
}

mean(resAfr["x",]); sd(resAfr["x",]) # x chr results
# -0.001864192 | 0.2250545
mean(resAfr["all",]); sd(resAfr["all",]) # all chrs results
# 0.003027238 | 0.2254082

# calculate empirical pvalue
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.2021212
sum(resAfr["all",]>=ca)/ncol(resAfr) # 0.2006; there are 1003 of 5000 greater
# test no random mating, ie assortative or disassortative
sum(abs(resAfr["all",])>=ca)/ncol(resAfr) # 0.39; there are 1950 of 5000 greater

(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.01230171
sum(resAfr["x",]>=cx)/ncol(resAfr) # 0.4244; there are 2122 of 5000 greater
# test no random mating, ie assortative or disassortative
sum(abs(resAfr["x",])>=cx)/ncol(resAfr) # 0.9582; there are 4791 of 5000 greater



###
## do for EUR
nms <- paste("chr",c(1:22,"X","All"),".local.CEU",sep="")
ethnData <- mates[,c(nms,"V2","sex")]
names(ethnData) <- c(1:22,"x","all","coriell.id","sex")

resEur <- getEst(ethnData,5432,withX=TRUE,allChrs=TRUE)

chrxF <- ethnData[ethnData$sex=="F",c(1:22,"x","all")]
chrxM <- ethnData[ethnData$sex=="M",c(1:22,"x","all")]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:24){
  if(i==23){filen <- paste("../assortative_mating/plots/localAncestry_Corr_allchr_mateAncest_CEU.pdf",sep="")
            mainTitle <- paste("Correlation of European local ancestry between mates\nacross all chromosomes")
            xlabAxis <- paste("Corr of All Chrs European local ancestry")
            cx <- corrChr["all"]
  }else if(i==24){filen <- paste("../assortative_mating/plots/localAncestry_Corr_xchr_mateAncest_CEU.pdf",sep="")
                  mainTitle <- paste("Correlation of European local ancestry between mates\non X chromosome")
                  xlabAxis <- paste("Corr of X Chr European local ancestry")    
                  cx <- corrChr["x"]
  }else{  
    filen <- paste("../assortative_mating/plots/localAncestry_Corr_",i,"chr_mateAncest_CEU.pdf",sep="")
    mainTitle <- paste("Correlation of European local ancestry between mates\non chromosome",i)
    xlabAxis <- paste("Corr of Chr",i,"European local ancestry")
    cx <- corrChr[i]
  }
  pdf(filen)
  hist(resEur[i,],main=mainTitle,xlab=xlabAxis,
       sub="Empirical distribution of 5000 resamples",breaks=50)
  abline(v=cx,col="red",lwd=2)
  cx <- format(cx,digits=3)
  legend("topright",(paste("Obs correlation =",cx)),col="red",
         lty=1,bg="white")
  dev.off()
}

mean(resEur["x",]); sd(resEur["x",]) # x chr results
#-0.000630983 | 0.2223752
mean(resEur["all",]); sd(resEur["all",]) # all chrs results
# 0.002695973 | 0.2249891

# calculate empirical pvalue
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.4498816
sum(resEur["all",]>=ca)/ncol(resEur) # 0.0196; there are 98 of 5000 greater
# test no random mating, ie assortative or disassortative
sum(abs(resEur["all",])>=ca)/ncol(resEur) # 0.0386; there are 193 of 5000 greater

(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.4796885
sum(resEur["x",]>=cx)/ncol(resEur) # 0.0164; there are 82 of 5000 greater
# test no random mating, ie assortative or disassortative
sum(abs(resEur["x",])>=cx)/ncol(resEur) # 0.029; there are 145 of 5000 greater


###
# save all results

save(res,file="../assortative_mating/local_HGDP_5000resamp.RData")
save(resEur,file="../assortative_mating/local_Euro_5000resamp.RData")
save(resAfr,file="../assortative_mating/local_YRI_5000resamp.RData")

#####
# are the males and females at the same proportion of ancestry?
mean(mates[mates$sex=="M","chrX.ceu"]-mates[mates$sex=="F","chrX.ceu"]) # -0.06013458
mean(mates[mates$sex=="M","chrX.hgdp"]-mates[mates$sex=="F","chrX.hgdp"]) # 0.07144417
mean(mates[mates$sex=="M","chrX.yri"]-mates[mates$sex=="F","chrX.yri"]) # -0.01130959

rm(list=ls())


#####
# 6. Plot genome wide local and ADMIX diff for all individs

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
ymin <- min(c(mxl$chrAll.local.CEU-mxl$chrAll.ceu,mxl$chrAll.local.YRI-mxl$chrAll.yri,
              mxl$chrAll.local.NAM-mxl$chrAll.hgdp))
ymax <- max(c(mxl$chrAll.local.CEU-mxl$chrAll.ceu,mxl$chrAll.local.YRI-mxl$chrAll.yri,
              mxl$chrAll.local.NAM-mxl$chrAll.hgdp))

pdf("../diff_localvsAdmix_autoWide.pdf",width=11)
plot(x=1:nrow(mxl),y=mxl$chrAll.local.CEU-mxl$chrAll.ceu,type="b",col="cyan",lwd=2,
     xlab="Individual",ylab="Average Local Ancestry - ADMIXTURE Estimate",main="Average Local Ancestry Autosomal-Wide
     Compared with ADMIXTURE Estimates Autosomal-Wide",ylim=c(ymin,ymax))
points(x=1:nrow(mxl),y=mxl$chrAll.local.YRI-mxl$chrAll.yri,type="b",col="magenta",lty=2,lwd=2)
points(x=1:nrow(mxl),y=mxl$chrAll.local.NAM-mxl$chrAll.hgdp,type="b",col="blue",lty=3,lwd=2)
legend("bottomleft",c(expression(paste("CEU, ",rho,"=0.993",sep="")),expression(paste("YRI, ",rho,"=0.893",sep="")),
                      expression(paste("NAm, ",rho,"=0.994",sep=""))),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)
abline(h=0,col="gray")
dev.off()

cor(mxl$chrAll.local.YRI,mxl$chrAll.yri) # 0.8932
cor(mxl$chrAll.local.CEU,mxl$chrAll.ceu) # 0.9935
cor(mxl$chrAll.local.NAM,mxl$chrAll.hgdp) # 0.9940

# make pairwise plots of the ancestries for each individ
pdf("../local_admix_YRI_pairs.pdf")
plot(mxl$chrAll.local.YRI,mxl$chrAll.yri,pch=19,col="magenta",xlab="Avg Autosomal Local Ancestry",
     ylab="Avg ADMIXTURE Estimate")
abline(a=0,b=1)
legend("topleft",c(expression(paste("YRI, ",rho,"=0.893",sep=""))),pch=19,col="magenta")
dev.off()

pdf("../local_admix_CEU_pairs.pdf")
plot(mxl$chrAll.local.CEU,mxl$chrAll.ceu,pch=19,col="cyan",xlab="Avg Autosomal Local Ancestry",
     ylab="Avg ADMIXTURE Estimate")
abline(a=0,b=1)
legend("topleft",c(expression(paste("CEU, ",rho,"=0.993",sep=""))),pch=19,col="cyan")
dev.off()

pdf("../local_admix_NAm_pairs.pdf")
plot(mxl$chrAll.local.NAM,mxl$chrAll.hgdp,pch=19,col="blue",xlab="Avg Autosomal Local Ancestry",
     ylab="Avg ADMIXTURE Estimate")
abline(a=0,b=1)
legend("topleft",c(expression(paste("NAm, ",rho,"=0.994",sep=""))),pch=19,col="blue")
dev.off()

rm(list=ls())


#####
# 7. Look at var of local ancestry est across individs

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))

varAcrossIndivids <- data.frame(matrix(NA,nrow=24,ncol=6))
colnames(varAcrossIndivids) <- c("ADMIX.yri","local.ancest.yri",
                                 "ADMIX.ceu","local.ancest.ceu",
                                 "ADMIX.nam","local.ancest.nam")

for(chr in 1:22){
  colAdmix <- paste("chr",chr,".yri",sep="")
  collocal <- paste("chr",chr,".local.YRI",sep="")
  varAcrossIndivids$ADMIX.yri[chr] <- var(mxl[,colAdmix])
  varAcrossIndivids$local.ancest.yri[chr] <- var(mxl[,collocal])
  
  colAdmix <- paste("chr",chr,".ceu",sep="")
  collocal <- paste("chr",chr,".local.CEU",sep="")
  varAcrossIndivids$ADMIX.ceu[chr] <- var(mxl[,colAdmix])
  varAcrossIndivids$local.ancest.ceu[chr] <- var(mxl[,collocal])
  
  colAdmix <- paste("chr",chr,".hgdp",sep="")
  collocal <- paste("chr",chr,".local.NAM",sep="")
  varAcrossIndivids$ADMIX.nam[chr] <- var(mxl[,colAdmix])
  varAcrossIndivids$local.ancest.nam[chr] <- var(mxl[,collocal])
}

# x chr now in row 23
  varAcrossIndivids$ADMIX.yri[23] <- var(mxl[,"chrX.yri"])
  varAcrossIndivids$local.ancest.yri[23] <- var(mxl[,"chrX.local.YRI"])

  varAcrossIndivids$ADMIX.ceu[23] <- var(mxl[,"chrX.ceu"])
  varAcrossIndivids$local.ancest.ceu[23] <- var(mxl[,"chrX.local.CEU"])
  
  varAcrossIndivids$ADMIX.nam[23] <- var(mxl[,"chrX.hgdp"])
  varAcrossIndivids$local.ancest.nam[23] <- var(mxl[,"chrX.local.NAM"])

# all chrs now in row 24
varAcrossIndivids$ADMIX.yri[24] <- var(mxl[,"chrAll.yri"])
varAcrossIndivids$local.ancest.yri[24] <- var(mxl[,"chrAll.local.YRI"])

varAcrossIndivids$ADMIX.ceu[24] <- var(mxl[,"chrAll.ceu"])
varAcrossIndivids$local.ancest.ceu[24] <- var(mxl[,"chrAll.local.CEU"])

varAcrossIndivids$ADMIX.nam[24] <- var(mxl[,"chrAll.hgdp"])
varAcrossIndivids$local.ancest.nam[24] <- var(mxl[,"chrAll.local.NAM"])

# looks like CEU has the largest difference in var looking over all the chrs

save(varAcrossIndivids,file="../assortative_mating/varAcrossIndivids_admixVsLocal.RData")
rm(list=ls())


#####
# 8. Plot genome wide local and ADMIX diff for all individs on X chr

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
ymin <- min(c(mxl$chrX.local.CEU-mxl$chrX.ceu,mxl$chrX.local.YRI-mxl$chrX.yri,
              mxl$chrX.local.NAM-mxl$chrX.hgdp))
ymax <- max(c(mxl$chrX.local.CEU-mxl$chrX.ceu,mxl$chrX.local.YRI-mxl$chrX.yri,
              mxl$chrX.local.NAM-mxl$chrX.hgdp))

cor(mxl$chrX.local.YRI,mxl$chrX.yri) # 0.9256
cor(mxl$chrX.local.CEU,mxl$chrX.ceu) # 0.9878
cor(mxl$chrX.local.NAM,mxl$chrX.hgdp) # 0.9898

pdf("../diff_localvsAdmix_Xchr.pdf",width=11)
plot(x=1:nrow(mxl),y=mxl$chrX.local.CEU-mxl$chrX.ceu,type="b",col="cyan",lwd=2,
     xlab="Individual",ylab="Average Local Ancestry - ADMIXTURE Estimate",main="Average Local Ancestry on X Chromosome
     Compared with ADMIXTURE Estimates on X Chromosome",ylim=c(ymin,ymax))
points(x=1:nrow(mxl),y=mxl$chrX.local.YRI-mxl$chrX.yri,type="b",col="magenta",lty=2,lwd=2)
points(x=1:nrow(mxl),y=mxl$chrX.local.NAM-mxl$chrX.hgdp,type="b",col="blue",lty=3,lwd=2)
legend("bottomleft",c(expression(paste("CEU, ",rho,"=0.988",sep="")),expression(paste("YRI, ",rho,"=0.926",sep="")),
                      expression(paste("NAm, ",rho,"=0.990",sep=""))),col=c("cyan","magenta","blue"),lty=c(1,2,3),lwd=2)
abline(h=0,col="gray")
dev.off()

# make pairwise plots of the ancestries for each individ
pdf("../local_admix_YRI_pairs_xchr.pdf")
plot(mxl$chrX.local.YRI,mxl$chrX.yri,pch=19,col="magenta",xlab="Avg X Chromosome Local Ancestry",
     ylab="ADMIXTURE X Chr Estimate")
abline(a=0,b=1)
legend("topleft",c(expression(paste("YRI, ",rho,"=0.9256",sep=""))),pch=19,col="magenta")
dev.off()

pdf("../local_admix_CEU_pairs_xchr.pdf")
plot(mxl$chrX.local.CEU,mxl$chrX.ceu,pch=19,col="cyan",xlab="Avg X Chromosome Local Ancestry",
     ylab="ADMIXTURE X Chr Estimate")
abline(a=0,b=1)
legend("topleft",c(expression(paste("CEU, ",rho,"=0.9878",sep=""))),pch=19,col="cyan")
dev.off()

pdf("../local_admix_NAm_pairs_xchr.pdf")
plot(mxl$chrX.local.NAM,mxl$chrX.hgdp,pch=19,col="blue",xlab="Avg X Chromosome Local Ancestry",
     ylab="ADMIXTURE X Chr Estimate")
abline(a=0,b=1)
legend("topleft",c(expression(paste("NAm, ",rho,"=0.9898",sep=""))),pch=19,col="blue")
dev.off()

rm(list=ls())


#####
# 9. Get mean (SD) for estimates by chromosome

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))

afr <- seq(from=80,to=146,by=3)
apply(mxl[,afr],2,function(x){print(paste(format(mean(x),digits=3)," (",format(sd(x),digits=3),")",sep=""))})

eur <- afr-1
apply(mxl[,eur],2,function(x){print(paste(format(mean(x),digits=3)," (",format(sd(x),digits=3),")",sep=""))})

nam <- afr+1
apply(mxl[,nam],2,function(x){print(paste(format(mean(x),digits=3)," (",format(sd(x),digits=3),")",sep=""))})


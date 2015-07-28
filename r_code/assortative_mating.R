#########################
# look at trios
# assess evidence of assortative mating

# function to resample 5000 times for each chr, x chr, and autosomes together
# calculates cor between randomly sampled mate pairs
getEst <- function(ethnData,chooseSeed){
  set.seed(chooseSeed)
  femList <- ethnData$coriell.id[ethnData$sex=="F"]
  malList <- ethnData$coriell.id[ethnData$sex=="M"]
  
  doOne <- function(){
    newF <- sample(femList,length(femList),replace=FALSE)
    newM <- sample(malList,length(malList),replace=FALSE)
    chrxF <- ethnData$x[match(newF,ethnData$coriell.id)]
    chrxM <- ethnData$x[match(newM,ethnData$coriell.id)]
    cx <- cor(chrxF,chrxM)
    
    newF <- sample(femList,length(femList),replace=FALSE)
    newM <- sample(malList,length(malList),replace=FALSE)
    chrAllF <- ethnData$all[match(newF,ethnData$coriell.id)]
    chrAllM <- ethnData$all[match(newM,ethnData$coriell.id)]
    ca <- cor(chrAllF,chrAllM)
    
    corrChr <- rep(NA,22)
    for(i in 1:22){
      newF <- sample(femList,length(femList),replace=FALSE)
      newM <- sample(malList,length(malList),replace=FALSE)
      chrF <- ethnData[match(newF,ethnData$coriell.id),i]
      chrM <- ethnData[match(newM,ethnData$coriell.id),i]
      corrChr[i] <- cor(chrF,chrM)
    }
    return(c(corrChr,ca,cx))
  }
  res <- replicate(5000,doOne())
  return(res)
}


mxl <- get(load("hapmap_mxlOnly_estimates.RData"))
table(table(mxl[,1]))
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
dim(mates) # 48 78

table(table(mates$superfamily))
#2  4  6 
#19  1  1

t <- table(mates$superfamily)
mates[is.element(mates$superfamily,names(t)[t==4]),]
# take NA19657, NA19658
mates <- mates[!is.element(mates[,2],c("NA19785","NA19786")),]
dim(mates)
mates[is.element(mates$superfamily,names(t)[t==6]),]
# take the two that have superfamilyUnrel==TRUE
mates <- mates[!is.element(mates[,2],c("NA19684","NA19661","NA19660","NA19685")),]
dim(mates) # 42 78; great!

table(table(mates[,1])) # 21 pairs
table(table(mates$superfamily)) # 21 pairs

mates <- merge(mates,parents,by.y="coriell.id",by.x="V2")
dim(mates); head(mates) # 42 84
# now these are all the mate pairs, w sex info
sum(duplicated(mates$coriell.id)) # 0, so all unique subjects

table(mates$family,mates$sex) # good, so one F & one M per family
mates <- mates[order(mates$family),]

## do for HGDP/Native Am
nms <- paste("chr",1:22,".hgdp",sep="")
ethnData <- mates[,c(nms,"V2","sex","chrX.hgdp","chrAll.hgdp")]
names(ethnData) <- c(1:22,"coriell.id","sex","x","all")

res <- getEst(ethnData,5432)
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.5364729
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.4672163
#(cx <- cor(chrxF,chrxM)) # 0.484148
#(ca <- cor(chrAllF,chrAllM)) # 0.484901

pdf("../assortative_mating/plots/Corr_Xchr_mateAncest.pdf")
hist(res[24,],main="Correlation of HGDP ancestry between mates\non the X chromosome",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of X chr HGDP ancestry",
     breaks=50)
abline(v=cx,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.536"))),col="red",
       lty=1,bg="white")
dev.off()

pdf("../assortative_mating/plots/Corr_allChr_mateAncest.pdf")
hist(res[23,],main="Correlation of HGDP ancestry between mates\nover all chromosomes",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of HGDP ancestry",
     breaks=100)
abline(v=ca,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.467"))),col="red",
       lty=1,bg="white")
dev.off()

chrxF <- ethnData[ethnData$sex=="F",1:22]
chrxM <- ethnData[ethnData$sex=="M",1:22]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:22){
  filen <- paste("../assortative_mating/plots/Corr_",i,"chr_mateAncest.pdf",sep="")
  mainTitle <- paste("Correlation of HGDP ancestry between mates\non chromosome",i)
  xlabAxis <- paste("Corr of Chr",i,"HGDP ancestry")
  cx <- corrChr[i]
  pdf(filen)
  hist(res[i,],main=mainTitle,xlab=xlabAxis,
       sub="Empirical distribution of 5000 resamples",breaks=50)
  abline(v=cx,col="red",lwd=2)
  cx <- format(cx,digits=3)
  legend("topright",(paste("Obs correlation =",cx)),col="red",
         lty=1,bg="white")
  dev.off()
}
mean(res[1,]); sd(res[1,]) # x chr results
#-0.00383267 | 0.2244598
mean(res[2,]); sd(res[2,]) # all chrs results
# -0.0002993807 |  0.220445

# calculate empirical pvalue
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.5364729
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.4672163

sum(res[1,]>=cx)/ncol(res) # 0.0068; there are 34 of 5000 greater
sum(res[2,]>=ca)/ncol(res) # 0.0172; there are 86 of 5000 greater

# get a ci based on the standard error=sqrt(p*(1-p)/5000
tmp <- sort(res[1,])
quantile(tmp,0.95) # 0.3651984 

tmp <- sort(res[2,])
quantile(tmp,0.95) # 0.3677846 

# calculate empirical pvalue for non-random mating
sum(abs(res[1,])>=cx)/ncol(res) # 0.0126; there are 63 of 5000 greater or smaller
sum(abs(res[2,])>=ca)/ncol(res) # 0.0322; there are 161 of 5000 greater or smaller


###
## do for YRI
nms <- paste("chr",1:22,".yri",sep="")
ethnData <- mates[,c(nms,"V2","sex","chrX.yri","chrAll.yri")]
names(ethnData) <- c(1:22,"coriell.id","sex","x","all")

resAfr <- getEst(ethnData,5432)
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.1471241
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.2611827

pdf("../assortative_mating/plots/Corr_Xchr_mateAncest_YRI.pdf")
hist(resAfr[24,],main="Correlation of African ancestry between mates\non the X chromosome",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of X chr YRI ancestry",
     breaks=50)
abline(v=cx,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.147"))),col="red",
       lty=1,bg="white")
dev.off()

pdf("../assortative_mating/plots/Corr_allChr_mateAncest_YRI.pdf")
hist(resAfr[23,],main="Correlation of African ancestry between mates\nover all chromosomes",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of YRI ancestry",
     breaks=50)
abline(v=ca,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.261"))),col="red",
       lty=1,bg="white")
dev.off()

chrxF <- ethnData[ethnData$sex=="F",1:22]
chrxM <- ethnData[ethnData$sex=="M",1:22]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:22){
  filen <- paste("../assortative_mating/plots/Corr_",i,"chr_mateAncest_YRI.pdf",sep="")
  mainTitle <- paste("Correlation of YRI ancestry between mates\non chromosome",i)
  xlabAxis <- paste("Corr of Chr",i,"YRI ancestry")
  cx <- corrChr[i]
  pdf(filen)
  hist(resAfr[i,],main=mainTitle,xlab=xlabAxis,
       sub="Empirical distribution of 5000 resamples",breaks=50)
  abline(v=cx,col="red",lwd=2)
  cx <- format(cx,digits=3)
  legend("topright",(paste("Obs correlation =",cx)),col="red",
         lty=1,bg="white")
  dev.off()
}

mean(resAfr[1,]); sd(resAfr[1,]) # x chr results
# -0.00167622 | 0.2235425
mean(resAfr[2,]); sd(resAfr[2,]) # all chrs results
# 0.004993024 | 0.2260868

# calculate empirical pvalue
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.1471241
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.2611827
sum(resAfr[1,]>=cx)/ncol(resAfr) # 0.2558; there are 1279 of 5000 greater
sum(resAfr[2,]>=ca)/ncol(resAfr) # 0.1392; there are 696 of 5000 greater

tmp <- sort(resAfr[1,])
quantile(tmp,0.95) # 0.3806093 

tmp <- sort(resAfr[2,])
quantile(tmp,0.95) # 0.385888  

# calculate empirical pvalue for non-random mating
sum(abs(resAfr[1,])>=cx)/ncol(resAfr) # 0.5302; there are 2651 of 5000 greater or smaller
sum(abs(resAfr[2,])>=ca)/ncol(resAfr) # 0.268; there are 1340 of 5000 greater or smaller



###
## do for EUR
nms <- paste("chr",1:22,".ceu",sep="")
ethnData <- mates[,c(nms,"V2","sex","chrX.ceu","chrAll.ceu")]
names(ethnData) <- c(1:22,"coriell.id","sex","x","all")

resEur <- getEst(ethnData,543287)
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.4855516
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.4781779

pdf("../assortative_mating/plots/Corr_Xchr_mateAncest_CEU.pdf")
hist(resEur[1,],main="Correlation of European ancestry between mates\non the X chromosome",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of X chr CEU ancestry",
     breaks=50)
abline(v=cx,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.486"))),col="red",
       lty=1,bg="white")
dev.off()

pdf("../assortative_mating/plots/Corr_allChr_mateAncest_CEU.pdf")
hist(resEur[2,],main="Correlation of European ancestry between mates\nover all chromosomes",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of CEU ancestry",
     breaks=50)
abline(v=ca,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.478"))),col="red",
       lty=1,bg="white")
dev.off()

chrxF <- ethnData[ethnData$sex=="F",1:22]
chrxM <- ethnData[ethnData$sex=="M",1:22]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:22){
  filen <- paste("../assortative_mating/plots/Corr_",i,"chr_mateAncest_CEU.pdf",sep="")
  mainTitle <- paste("Correlation of CEU ancestry between mates\non chromosome",i)
  xlabAxis <- paste("Corr of Chr",i,"CEU ancestry")
  cx <- corrChr[i]
  pdf(filen)
  hist(resEur[i,],main=mainTitle,xlab=xlabAxis,
       sub="Empirical distribution of 5000 resamples",breaks=50)
  abline(v=cx,col="red",lwd=2)
  cx <- format(cx,digits=3)
  legend("topright",(paste("Obs correlation =",cx)),col="red",
         lty=1,bg="white")
  dev.off()
}

mean(resEur[1,]); sd(resEur[1,]) # x chr results
# -0.001262082 | 0.2223165
mean(resEur[2,]); sd(resEur[2,]) # all chrs results
# 3.818851e-05 | 0.2217552

# calculate empirical pvalue
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # 0.4855516
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.4781779
sum(resEur[1,]>=cx)/ncol(resEur) # 0.0114; there are 57 of 5000 greater
sum(resEur[2,]>=ca)/ncol(resEur) # 0.0146; there are 73 of 5000 greater

tmp <- sort(resEur[1,])
quantile(tmp,0.95) # 0.3495509 

tmp <- sort(resEur[2,])
quantile(tmp,0.95) # 0.3403861 

# calculate empirical pvalue for non-random mating
sum(abs(resEur[1,])>=cx)/ncol(resEur) # 0.024; there are 120 of 5000 greater or smaller
sum(abs(resEur[2,])>=ca)/ncol(resEur) # 0.0278; there are 139 of 5000 greater or smaller


###
# save all results

save(res,file="../assortative_mating/HGDP_5000resamp.RData")
save(resEur,file="../assortative_mating/Euro_5000resamp.RData")
save(resAfr,file="../assortative_mating/YRI_5000resamp.RData")

#####
# are the males and females at the same proportion of ancestry?
mean(mates[mates$sex=="M","chrX.ceu"]-mates[mates$sex=="F","chrX.ceu"]) # -0.06013458
mean(mates[mates$sex=="M","chrX.hgdp"]-mates[mates$sex=="F","chrX.hgdp"]) # 0.07144417
mean(mates[mates$sex=="M","chrX.yri"]-mates[mates$sex=="F","chrX.yri"]) # -0.01130959


rm(list=ls())

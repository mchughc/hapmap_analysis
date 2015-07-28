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

asw <- get(load("hapmap_aswOnly_estimates.RData"))
table(table(asw[,1]))
#1  2  3 
#9 21 12 

# merge in sex info
ped <- read.csv("~/Dropbox/asw_pedigree.csv",header=TRUE,as.is=TRUE)
all(is.element(asw[,2],ped$individ)) # TRUE
ped <- ped[is.element(ped$individ,asw[,2]),]
dim(ped) # 87 10
table(table(ped$family)) # only 11 trios. dang.
t <- table(ped$family)
trs <- names(t)[t==3]
parents <- ped[is.element(ped$family,trs)&ped$mother==0&ped$father==0,]
dim(parents); head(parents) # 22 10; good
table(table(parents$family)) # 11 of 2

mates <- asw[is.element(asw$V2,parents$individ),]
dim(mates) # 22 76

table(table(mates$V1))
#2 
#11

mates <- merge(mates,parents,by.y="individ",by.x="V2")
dim(mates); head(mates) # 22 85
# now these are all the mate pairs, w sex info
sum(duplicated(mates$V2)) # 0, so all unique subjects

table(mates$family,mates$sex) # good, so one F & one M per family
mates <- mates[order(mates$family),]

## do for HGDP/Native Am
nms <- paste("chr",1:22,".hgdp",sep="")
ethnData <- mates[,c(nms,"V2","sex","chrX.hgdp","chrAll.hgdp")]
names(ethnData) <- c(1:22,"coriell.id","sex","x","all")

res <- getEst(ethnData,5432)
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # -0.05285695
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.2144943
#(cx <- cor(chrxF,chrxM)) # 0.484148
#(ca <- cor(chrAllF,chrAllM)) # 0.484901

pdf("../assortative_mating/plots/asw_Corr_Xchr_mateAncest.pdf")
hist(res[24,],main="Correlation of HGDP ancestry between mates\non the X chromosome",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of X chr HGDP ancestry",
     breaks=50)
abline(v=cx,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=-0.053"))),col="red",
       lty=1,bg="white")
dev.off()

pdf("../assortative_mating/plots/asw_Corr_allChr_mateAncest.pdf")
hist(res[23,],main="Correlation of HGDP ancestry between mates\nover all chromosomes",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of HGDP ancestry",
     breaks=100)
abline(v=ca,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.214"))),col="red",
       lty=1,bg="white")
dev.off()

chrxF <- ethnData[ethnData$sex=="F",1:22]
chrxM <- ethnData[ethnData$sex=="M",1:22]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:22){
  filen <- paste("../assortative_mating/plots/asw_Corr_",i,"chr_mateAncest.pdf",sep="")
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
# 0.001104364 | 0.3144921
mean(res[2,]); sd(res[2,]) # all chrs results
# -0.002637127 | 0.3140781

# calculate empirical pvalue
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # -0.05285695
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.2144943

sum(res[1,]>=cx)/ncol(res) # 0.5636; there are 2818 of 5000 greater
sum(res[2,]>=ca)/ncol(res) # 0.2344; there are 1172 of 5000 greater

# get a ci based on the standard error=sqrt(p*(1-p)/5000
tmp <- sort(res[1,])
quantile(tmp,0.95) # 0.5345474 

tmp <- sort(res[2,])
quantile(tmp,0.95) # 0.6039933 

# calculate empirical pvalue for non-random mating
sum(abs(res[1,])>=cx)/ncol(res) # 1; there are 5000 of 5000 greater or smaller
sum(abs(res[2,])>=ca)/ncol(res) # 0.5316; there are 2658 of 5000 greater or smaller


###
## do for YRI
nms <- paste("chr",1:22,".yri",sep="")
ethnData <- mates[,c(nms,"V2","sex","chrX.yri","chrAll.yri")]
names(ethnData) <- c(1:22,"coriell.id","sex","x","all")

resAfr <- getEst(ethnData,5432)
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # -0.3348569
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.05435167

pdf("../assortative_mating/plots/asw_Corr_Xchr_mateAncest_YRI.pdf")
hist(resAfr[24,],main="Correlation of African ancestry between mates\non the X chromosome",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of X chr YRI ancestry",
     breaks=50)
abline(v=cx,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=-0.335"))),col="red",
       lty=1,bg="white")
dev.off()

pdf("../assortative_mating/plots/asw_Corr_allChr_mateAncest_YRI.pdf")
hist(resAfr[23,],main="Correlation of African ancestry between mates\nover all chromosomes",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of YRI ancestry",
     breaks=50)
abline(v=ca,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.054"))),col="red",
       lty=1,bg="white")
dev.off()

chrxF <- ethnData[ethnData$sex=="F",1:22]
chrxM <- ethnData[ethnData$sex=="M",1:22]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:22){
  filen <- paste("../assortative_mating/plots/asw_Corr_",i,"chr_mateAncest_YRI.pdf",sep="")
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
# -0.002484399 | 0.3172723
mean(resAfr[2,]); sd(resAfr[2,]) # all chrs results
# -0.0001806461 | 0.3191402

# calculate empirical pvalue
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # -0.3348569
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.05435167
sum(resAfr[1,]>=cx)/ncol(resAfr) # 0.8424; there are 4212 of 5000 greater
sum(resAfr[2,]>=ca)/ncol(resAfr) # 0.3648; there are 1824 of 5000 greater

tmp <- sort(resAfr[1,])
quantile(tmp,0.95) # 0.5455219

tmp <- sort(resAfr[2,])
quantile(tmp,0.95) # 0.62289   

# calculate empirical pvalue for non-random mating
sum(abs(resAfr[1,])>=cx)/ncol(resAfr) # 1; there are 5000 of 5000 greater or smaller
sum(abs(resAfr[2,])>=ca)/ncol(resAfr) # 0.8712; there are 4356 of 5000 greater or smaller



###
## do for EUR
nms <- paste("chr",1:22,".ceu",sep="")
ethnData <- mates[,c(nms,"V2","sex","chrX.ceu","chrAll.ceu")]
names(ethnData) <- c(1:22,"coriell.id","sex","x","all")

resEur <- getEst(ethnData,543287)
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # -0.2718763
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.05018522

pdf("../assortative_mating/plots/asw_Corr_Xchr_mateAncest_CEU.pdf")
hist(resEur[1,],main="Correlation of European ancestry between mates\non the X chromosome",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of X chr CEU ancestry",
     breaks=50)
abline(v=cx,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=-0.272"))),col="red",
       lty=1,bg="white")
dev.off()

pdf("../assortative_mating/plots/asw_Corr_allChr_mateAncest_CEU.pdf")
hist(resEur[2,],main="Correlation of European ancestry between mates\nover all chromosomes",
     sub="Empirical distribution of 5000 resamples",xlab="Corr of CEU ancestry",
     breaks=50)
abline(v=ca,col="red",lwd=2)
legend("topright",c(expression(paste("Obs correlation ",rho,"=0.050"))),col="red",
       lty=1,bg="white")
dev.off()

chrxF <- ethnData[ethnData$sex=="F",1:22]
chrxM <- ethnData[ethnData$sex=="M",1:22]
corrChr <- diag(cor(chrxF,chrxM))
for(i in 1:22){
  filen <- paste("../assortative_mating/plots/asw_Corr_",i,"chr_mateAncest_CEU.pdf",sep="")
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
# -0.008317919 | 0.3130222
mean(resEur[2,]); sd(resEur[2,]) # all chrs results
# 0.003013694 | 0.3197325

# calculate empirical pvalue
(cx <- cor(ethnData$x[ethnData$sex=="F"],ethnData$x[ethnData$sex=="M"])) # -0.2718763
(ca <- cor(ethnData$all[ethnData$sex=="F"],ethnData$all[ethnData$sex=="M"])) # 0.05018522
sum(resEur[1,]>=cx)/ncol(resEur) # 0.7876; there are 3938 of 5000 greater
sum(resEur[2,]>=ca)/ncol(resEur) # 0.3876; there are 1938 of 5000 greater

tmp <- sort(resEur[1,])
quantile(tmp,0.95) # 0.5448346  

tmp <- sort(resEur[2,])
quantile(tmp,0.95) # 0.6137011  

# calculate empirical pvalue for non-random mating
sum(abs(resEur[1,])>=cx)/ncol(resEur) # 1; there are 5000 of 5000 greater or smaller
sum(abs(resEur[2,])>=ca)/ncol(resEur) # 0.8884; there are 4442 of 5000 greater or smaller


###
# save all results

save(res,file="../assortative_mating/asw_HGDP_5000resamp.RData")
save(resEur,file="../assortative_mating/asw_Euro_5000resamp.RData")
save(resAfr,file="../assortative_mating/asw_YRI_5000resamp.RData")

#####
# are the males and females at the same proportion of ancestry?
mean(mates[mates$sex=="M","chrX.ceu"]-mates[mates$sex=="F","chrX.ceu"]) # -0.03274103
mean(mates[mates$sex=="M","chrX.hgdp"]-mates[mates$sex=="F","chrX.hgdp"]) # -0.008744491
mean(mates[mates$sex=="M","chrX.yri"]-mates[mates$sex=="F","chrX.yri"]) # 0.04148545

rm(list=ls())

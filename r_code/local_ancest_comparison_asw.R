### compare the ADMIXTURE results for each individual, by chromosome, to the average of the local ancestry
## for the HapMap ASW samples

## Contents:
# 1. Process files
# 2. Redo assortative mating analysis with local admixture estimates
# 3. Get mean (SD) for estimates by chromosome




#####
# 1. Process files

asw <- get(load("hapmap_aswOnly_estimates.RData"))

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


# get average local ancestry for MXL individuals by chromosome
for(i in 1:23){
  ancest <- read.table(paste("../assortative_mating/ASW_ancestry/HAPMAP_ASW_CHR",i,".txt",sep=""),header=T)
  stopifnot(colnames(ancest)==c("FID","IID","CEU","YRI","NAM"))
  colnames(ancest) <- c("FID","IID",paste("chr",i,".local.CEU",sep=""),paste("chr",i,".local.YRI",sep=""),
                        paste("chr",i,".local.NAM",sep=""))
  asw <- merge(asw,ancest[,2:5],by.x="V2",by.y="IID",all.x=TRUE)
}

colnames(asw)[colnames(asw)=="chr23.local.CEU"] <- "chrX.local.CEU"
colnames(asw)[colnames(asw)=="chr23.local.YRI"] <- "chrX.local.YRI"
colnames(asw)[colnames(asw)=="chr23.local.NAM"] <- "chrX.local.NAM"

# get average autosomal estimates too
ceuCols <- seq(from=77,to=142,by=3)
colnames(asw)[ceuCols]
asw$chrAll.local.CEU <- rowMeans(asw[,c(ceuCols)])

yriCols <- ceuCols+1
colnames(asw)[yriCols]
asw$chrAll.local.YRI <- rowMeans(asw[,c(yriCols)])

namCols <- yriCols+1
colnames(asw)[namCols]
asw$chrAll.local.NAM <- rowMeans(asw[,c(namCols)])

save(asw,file="asw_estimates_withLocal.RData")

rm(list=ls())


#####
# 2. Redo assortative mating analysis with local admixture estimates

source("getEst.R")

asw <- get(load("asw_estimates_withLocal.RData"))
table(table(asw[,2]))
#  1  2  3 
#  9 21 12 

# merge in sex info
ped <- read.csv("~/Dropbox/asw_pedigree.csv",header=TRUE,as.is=TRUE)
all(is.element(asw$V2,ped$individ)) # TRUE
ped <- ped[is.element(ped$individ,asw$V2),]
dim(ped) # 87 10
table(table(ped$family)) # only 11 trios. dang.
t <- table(ped$family)
trs <- names(t)[t==3]
parents <- ped[is.element(ped$family,trs)&ped$mother==0&ped$father==0,]
dim(parents); head(parents) # 22 10; good
table(table(parents$family)) # 11 of 2

mates <- asw[is.element(asw$V2,parents$individ),]
dim(mates) # 22 148
table(table(mates$V1))
#2 
#11

mates <- merge(mates,parents,by.y="individ",by.x="V2")
dim(mates); head(mates) # 22 157
# now these are all the mate pairs, w sex info
sum(duplicated(mates$V2)) # 0, so all unique subjects

table(mates$family,mates$sex) # good, so one F & one M per family
mates <- mates[order(mates$family),]

###
## do for HGDP/Native Am
nms <- paste("chr",c(1:22,"X","All"),".local.NAM",sep="")
ethnData <- mates[,c(nms,"V2","sex")]
names(ethnData) <- c(1:22,"x","all","coriell.id","sex")

# hmm, a lot of zeroes
dim(ethnData) # 22 26; so a 24 would mean all entries equal zero
apply(ethnData[,1:24],1,function(x){sum(x==0)})
# no samples have identically zero for all chrs in NAm ancestry estimation
# one sample has 22 0's though -- sample 15?
ethnData[15,] # chr 18 is only non-zero, and avg is thus slightly larger than 0

# don't go any further on these analyses

rm(list=ls())


#####
# 3. Get mean (SD) for estimates by chromosome

asw <- get(load("asw_estimates_withLocal.RData"))

afr <- seq(from=78,to=144,by=3)
apply(asw[,afr],2,function(x){print(paste(format(mean(x),digits=3)," (",format(sd(x),digits=3),")",sep=""))})

eur <- afr-1
apply(asw[,eur],2,function(x){print(paste(format(mean(x),digits=3)," (",format(sd(x),digits=3),")",sep=""))})

nam <- afr+1
apply(asw[,nam],2,function(x){print(paste(format(mean(x),digits=3)," (",format(sd(x),digits=3),")",sep=""))})

rm(list=ls())


#####
# 4. Get diff of X chr ancestry between mates

males <- mates[mates$sex=="M",]
females <- mates[mates$sex=="F",]

stopifnot(males$gcc.family==females$gcc.family)

range(males$chrX.local.NAM-females$chrX.local.NAM)
range(males$chrX.local.YRI-females$chrX.local.YRI)
range(males$chrX.local.CEU-females$chrX.local.CEU)

# combine all proportions for each chr into one file
chr1 <- read.table("~/Dropbox/Tim_Caitlin_X/HapMap_ASW_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_ASW_CHR1_result.txt",
                   header=F,as.is=T)
dim(chr1) # 518 6
chr1 <- data.frame(chr1)
chr2 <- read.table("~/Dropbox/Tim_Caitlin_X/HapMap_ASW_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_ASW_CHR2_result.txt",
                   header=F,as.is=T)
dim(chr2) # 518 6
chr2 <- data.frame(chr2)

all(chr1[,1]==chr2[,1]) # TRUE
all(chr1[,2]==chr2[,2]) # TRUE
allchrs <- chr1
for(i in 2:22){
  filen <- paste("~/Dropbox/Tim_Caitlin_X/HapMap_ASW_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_ASW_CHR",i,"_result.txt",sep="")
  newc <- read.table(filen,header=F,as.is=T)
  stopifnot(allchrs[,1]==newc[,1])
  stopifnot(allchrs[,2]==newc[,2])
  names(newc)[4:6] <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))
  allchrs <- cbind(allchrs,newc[,4:6])
}
filen <- "~/Dropbox/Tim_Caitlin_X/HapMap_ASW_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_CHRX_ONE_MALE_ALLELE_0_ASW_result.txt"
newc <- read.table(filen,header=F,as.is=T)
stopifnot(allchrs[,1]==newc[,1])
stopifnot(allchrs[,2]==newc[,2])
names(newc)[4:6] <- c("chrX.ceu","chrX.yri","chrX.hgdp")
allchrs <- cbind(allchrs,newc[,4:6])
names(allchrs)[4:6] <- c("chr1.ceu","chr1.yri","chr1.hgdp")

# no all autosomal SNPs analysis for the ASW samples
# filen <- "~/Dropbox/Tim_Caitlin_X/HapMap_ASW_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_USING_ALL_AUTOSOMAL_CHROMOSOMES_result.txt"
# newc <- read.table(filen,header=F,as.is=T)
# stopifnot(allchrs[,1]==newc[,1])
# stopifnot(allchrs[,2]==newc[,2])
# names(newc)[4:6] <- c("chrAll.ceu","chrAll.yri","chrAll.hgdp")
# allchrs <- cbind(allchrs,newc[,4:6])

ceuCols <- seq(from=4,to=69,by=3)
names(allchrs)[ceuCols]; names(allchrs)[ceuCols+1]; names(allchrs)[ceuCols+2]
allchrs$chrAll.ceu <- apply(allchrs[,ceuCols],1,mean)
allchrs$chrAll.yri <- apply(allchrs[,(ceuCols+1)],1,mean)
allchrs$chrAll.hgdp <- apply(allchrs[,(ceuCols+2)],1,mean)

table(allchrs[,4]) # so there are 266 hgdp subj, 165 yri subj and 87 asw subj
# only 45 are unrelated
asw <- allchrs[!is.element(allchrs[,4],c(0,1)),]
dim(asw) # 87 75

# add an unrelated variable
fam <- read.table("~/Dropbox/Tim_Caitlin_X/HapMap_ASW_Ancestry_Estimates/HapMap_ASW_Unrelated_Set_of_45.output")
sum(is.element(fam[,2],asw[,2])) # 45
asw$unrelated <- is.element(asw[,2],fam[,2])
table(asw$unrelated,exclude=NULL) # 45 true; 42 false

save(asw,file="hapmap_aswOnly_estimates.RData")
save(allchrs,file="hapmap_asw_estimates_allchrs.RData")

rm(list=ls())

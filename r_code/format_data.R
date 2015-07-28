# combine all proportions for each chr into one file
chr1 <- read.table("~/Dropbox/HapMap_MXL_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_CHR1_result.txt",
                   header=F,as.is=T)
dim(chr1) # 517 6
chr1 <- data.frame(chr1)
chr2 <- read.table("~/Dropbox/HapMap_MXL_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_CHR2_result.txt",
                   header=F,as.is=T)
dim(chr2) # 517 6
chr2 <- data.frame(chr2)

all(chr1[,1]==chr2[,1]) # TRUE
all(chr1[,2]==chr2[,2]) # TRUE
allchrs <- chr1
for(i in 2:22){
  filen <- paste("~/Dropbox/HapMap_MXL_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_CHR",i,"_result.txt",sep="")
  newc <- read.table(filen,header=F,as.is=T)
  stopifnot(allchrs[,1]==newc[,1])
  stopifnot(allchrs[,2]==newc[,2])
  names(newc)[4:6] <- c(paste("chr",i,".ceu",sep=""),paste("chr",i,".yri",sep=""),paste("chr",i,".hgdp",sep=""))
  allchrs <- cbind(allchrs,newc[,4:6])
}
filen <- "~/Dropbox/HapMap_MXL_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_CHRX_ONE_MALE_ALLELE_0_result.txt"
newc <- read.table(filen,header=F,as.is=T)
stopifnot(allchrs[,1]==newc[,1])
stopifnot(allchrs[,2]==newc[,2])
names(newc)[4:6] <- c("chrX.ceu","chrX.yri","chrX.hgdp")
allchrs <- cbind(allchrs,newc[,4:6])
names(allchrs)[4:6] <- c("chr1.ceu","chr1.yri","chr1.hgdp")

filen <- "~/Dropbox/HapMap_MXL_Ancestry_Estimates/FRAPPE_MERGED_HGDP_HAPMAP_RELEASE3_PLINK_FILE_USING_ALL_AUTOSOMAL_CHROMOSOMES_result.txt"
newc <- read.table(filen,header=F,as.is=T)
stopifnot(allchrs[,1]==newc[,1])
stopifnot(allchrs[,2]==newc[,2])
names(newc)[4:6] <- c("chrAll.ceu","chrAll.yri","chrAll.hgdp")
allchrs <- cbind(allchrs,newc[,4:6])

table(allchrs[,4]) # so there are 266 hgdp subj, 165 yri subj and 86 mxl subj
# only 53 are unrelated
mxl <- allchrs[!is.element(allchrs[,4],c(0,1)),]
dim(mxl) # 86 75

# add an unrelated variable
fam <- read.table("~/Dropbox/HapMap_MXL_Ancestry_Estimates/HapMap_MXL_Unrelated_Set_of_53.output")
sum(is.element(fam[,2],mxl[,2])) # 53
mxl$unrelated <- is.element(mxl[,2],fam[,2])
table(mxl$unrelated,exclude=NULL) # 53 true; 33 false

## does this include found cryptic relatedness from tim's paper?
# add a new col of 'superfamily'
mxl$superfamily <- mxl[,1] # start as the family id
mxl$superfamily[is.element(mxl[,1],c("M012","M008","M011","2382"))] <- "M012"
mxl$superfamily[is.element(mxl[,1],c("M032","M007"))] <- "M032"
length(unique(mxl$superfamily)) # 29
length(unique(mxl[,1])) # 33
# ok, good

# update the unrelated variable
mxl$superfamUnrel <- mxl$unrelated
mxl[is.element(mxl$superfamily,"M012"),]
mxl$superfamUnrel[is.element(mxl$superfamily,"M012")] <- FALSE
mxl$superfamUnrel[is.element(mxl[,2],c("NA19663","NA19664"))] <- TRUE

mxl[is.element(mxl$superfamily,"M032"),]
# ok, this is ok - the 3rd unrel sample is related at 4th degree, so it's ok

save(mxl,file="hapmap_mxlOnly_estimates.RData")
save(allchrs,file="hapmap_mxl_estimates_allchrs.RData")

rm(list=ls())


setwd("/projects/geneva/geneva_sata/HapMap/HapMap3_r3/sample_snp_annot")
b36 <- get(load("hapmap3_r3_b36.1Msnp.v4.RData"))
table(b36$chromosome) # 13273 x chr snps

b37 <- get(load("hapmap3_r3_b37.1Msnp.v6.RData"))
dim(b37) # 936978 10
table(b37$chromosome) # dang! only 2 chr 23 snps????

# we have recomb rates for x chr but it's in build 37
# need to backtrack to build 36

recomb37 <- read.table("/projects/geneva/geneva_sata/caitlin/genetic_map_GRCh37_chrX.txt",header=TRUE,as.is=TRUE)
head(recomb37); dim(recomb37) # 89589 4

# need to merge into the recomb37 file the build 36 chr and position
names(recomb37) <- c("Chromosome.37","Position.bp.37","Rate.cM.Mb.37","Map.cM.37")

# read in the list of rsIDs from build 36 that lisa sent me
rsb36 <- read.table("/projects/geneva/geneva_sata/caitlin/X_CHR_rs_IDs_b36_HAPMAP_CEU_YRI_MEX_and_HGDP_NATIVE_AMERICANS.txt",
                    header=TRUE,as.is=TRUE)
dim(rsb36) # 5685 1

sum(is.element(rsb36[,1],b36$rs.id)) # 5594
rsb36 <- merge(rsb36,pData(b36)[,c("chromosome","rs.id","position")],by.x="rsID",by.y="rs.id",all.x=TRUE)
dim(rsb36); head(rsb36) # 5685 3 | 3 cols, rsID, chr, position

rm(list=ls())


#####
# take build 37 recomb rate file and liftover to build 36, then merge with lisa's list of rsIDs and b36 positions

library(rtracklayer)

# read in recomb rates for build 37; de-update that to build 36
recomb37 <- read.table("/projects/geneva/geneva_sata/caitlin/genetic_map_GRCh37_chrX.txt",header=TRUE,as.is=TRUE)
head(recomb37); dim(recomb37) # 89589 4

snpAnnot <- SnpAnnotationDataFrame(data=data.frame(rsID=paste("rs",1:nrow(recomb37),sep=""),
                                  snpID=1:nrow(recomb37),position=recomb37$Position.bp.,
                                  chromosome=rep(as.integer(23),nrow(recomb37))))

# do liftover to build 36, get new build 36 positions
# send those back to lisa and see if that's sufficient
chain.file <- "/projects/geneva/gcc-fs2/SNP_annotation/UCSC_downloads/hg19/liftOver/hg19ToHg18.over.chain"
converted <- convertBuild(snpAnnot,chain.file)
# now converted has build 36 positions

# save the build 36 positions
lifts <- converted[,c("position","position.converted")]
colnames(lifts) <- c("b37.positions","b36.position")

recomb37 <- merge(recomb37,lifts,by.x="Position.bp.",by.y="b37.positions")
colnames(recomb37)[1] <- "b37.position"

write.table(recomb37,file="/projects/geneva/geneva_sata/caitlin/b37_b36_recomb_xchr.txt",quote=FALSE,col.names=TRUE)
# send this file to lisa
rm(list=ls())


#####








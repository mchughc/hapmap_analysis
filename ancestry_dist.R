dat <- read.table("WHI_AA_X_AUTOSOMAL_ANCESTRY/HUA_WHI_ANCESTRY_DISTRIBUTION_ALL_08_2010.txt",header=T,as.is=T,sep="\t")
dim(dat) # 12008 5

# so there are 12008 subjects
# want the proportion of african per subject, across the autosomes, and for the xchr

xinfo <- read.table("WHI_AA_X_AUTOSOMAL_ANCESTRY/sampleinfoX.txt",header=T,as.is=T)
# 8421 4
table(xinfo$dropMe,exclude=NULL)
#FALSE  TRUE  <NA> 
#  8150   271     0 

hist(as.numeric(xinfo$AfrIX[xinfo$dropMe==FALSE]))
hist(dat$YRI)
sum(is.element(dat$FID,xinfo$famID)) # 8421
all(dat$FID==xinfo$famID) # TRUE

dat <- dat[is.element(dat$FID,xinfo$famID),]
dim(dat) # 8421 5

ancest <- cbind(xinfo,dat[,"YRI"])
head(ancest); dim(ancest) # 8421 5
names(ancest) <- c("ID","famID","X","dropMe","Auto")

ancestf <- ancest[!ancest$dropMe,]
dim(ancestf) # 8150 5
ancestf$XAutoDiff <- ancestf$Auto-ancestf$X
ancestf$diff_ind <- ancestf$XAutoDiff>0

var(ancestf$X) # 0.03342007
var(ancestf$Auto) # 0.01910483

### testing
t.test(ancestf$X,ancestf$Auto)
#t = 18.0128, df = 15171.1, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0 
#95 percent confidence interval:
#  0.04075228 0.05070442 
#sample estimates:
#  mean of x mean of y 
#0.8128387 0.7671104

sum(ancestf$diff_ind)/nrow(ancestf) # 0.3193
# under the null, this proportion should be 1/2

### data summaries
var(ancestf$XAutoDiff) # 0.01977044
mean(ancestf$XAutoDiff) # -0.04572835

png("hist_autoXdiff.png")
hist(ancestf$XAutoDiff,xlab="Diff between Autosomal & X Chr",main="Difference between African Ancestry on Autosomes and X Chr")
dev.off()

png("xchr_AutoProp.png")
plot(ancestf$Auto,ancestf$X,xlab="Autosomal",ylab="X chr",sub="line is X=Y",
     main="Proportion of Afr Ancestry in 8150 Subjects",col="#FF000060",xlim=c(0,1),ylim=c(0,1))
abline(0,1,lty=2)
dev.off()


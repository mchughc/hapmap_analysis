#####
## get empirical distribution

setwd("~/Documents/tim_stuff/r_code")
source("emp_fisher.R")
library(CAnD)

mxl <- get(load("mxl_estimates_withLocal_v2.RData"))
mxl <- mxl[mxl$unrelated,]
dim(mxl) # 53 150

niters <- 1e6
ceu <- emp_fisher(mxl[,seq(from=79,to=146,by=3)],niters)
nam <- emp_fisher(mxl[,seq(from=81,to=147,by=3)],niters)
yri <- emp_fisher(mxl[,seq(from=80,to=147,by=3)],niters)

# get the actual estimates for the dat
true_ceu <- CAnD(mxl[,seq(from=79,to=146,by=3)],bonfCorr=FALSE)
true_nam <- CAnD(mxl[,seq(from=81,to=147,by=3)],bonfCorr=FALSE)
true_yri <- CAnD(mxl[,seq(from=80,to=147,by=3)],bonfCorr=FALSE)

summary(ceu)
summary(nam)
summary(yri)

write.table(ceu,file="ceu_1Miters.txt",quote=FALSE,row.names=FALSE)
write.table(nam,file="nam_1Miters.txt",quote=FALSE,row.names=FALSE)
write.table(yri,file="yri_1Miters.txt",quote=FALSE,row.names=FALSE)

q("no")


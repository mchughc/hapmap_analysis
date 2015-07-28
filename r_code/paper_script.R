
## script to process simulation results

# make a plot of the type I error results?

param <- c(0.0118,0.0053,0.0004,)
res <- read.csv("Dropbox/tim_stuff/paired_ttest_res_allRuns_01.csv",header=TRUE,as.is=TRUE)
dim(res) # 500 7

Parametric CAnD test	Non-parametric CAnD test
0.01	0.0118	0.0066
0.005	0.0053	0.0025
0.001	0.0004	0.0012
0.0001	0.0000	0.0000


# make a plot of the table 3 pvalues, for mxl and asw independently
# 3 dots at each ancestry for all samples, just males, and just females
# xaxis=afr, euro, Nam
# yaxis=-log10(pvalue)

res <- data.frame(Ancestry=rep(c("African","European","Native American"),times=3),set=c(rep("All Samples",3),rep("Females",3),rep("Males",3)))
res <- rbind(res,res)
res$HapMap <- c(rep("MXL",9),rep("ASW",9))
res$pvalue <- c(0.027,8.171e-05,2.248e-05,0.248,0.006,0.001,0.076,0.009,0.009,
                0.135,0.072,0.036,0.690,0.690,0.015,0.115,0.041,0.823)

library(ggplot2)
pdf("Documents/tim_stuff/ancestry_differences/plots/non_param_cand_pvalues.pdf")
ggplot(res,aes(x=Ancestry,y=-log10(pvalue),color=set,shape=HapMap)) + geom_point()
dev.off()


## do with local ancestry now
res <- data.frame(Ancestry=rep(c("African","European","Native American"),times=3),set=c(rep("All Samples",3),rep("Females",3),rep("Males",3)))
res <- rbind(res,res)
head(res)
tail(res)
res$HapMap <- c(rep("MXL",9),rep("ASW",9))
head(res)
res$pvalue <- c(0.05343881,5.550521e-06, 5.550521e-06, # afr, euro, nam in all MXL samp
0.4420683, 0.00151372, 0.00151372, # afr, euro, nam in female MXL samp
0.0755187, 0.002493918, 0.002493918, # afr euro nam in male MXL samp
0.2326932,0.1351565,0.0008240824, # afr euro nam in all ASW samp
1.0000001,0.690038,0.1077521, # afr euro nam in female ASW samp
0.1153183, 0.1153183, 0.002576828# afr euro nam in male ASW samp
)
head(res)
res
library(ggplot2)
pdf("../ancestry_differences/plots/non_param_cand_pvalues_localAncest.pdf")
ggplot(res,aes(x=Ancestry,y=-log10(pvalue),color=set,shape=HapMap)) + geom_point() + theme_bw()
dev.off()


#### 
# just plot the MXL
pdf("../ancestry_differences/plots/mxl_non_param_cand_pvalues_localAncest.pdf")
ggplot(res[res$HapMap=="MXL",],aes(x=Ancestry,y=-log10(pvalue),color=set, shape=set)) + geom_point(size=3) + theme_bw()
dev.off()


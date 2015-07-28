
# function to create an empirical distribution of fisher's z-scores
# resample sets of 23 chrs, do cand and calculate p-value
# save only pvalue

# this function does niters
emp_fisher <- function(ancestProp,niters){
  # ancestProp is a df with nrow=nsample x ncol=23 
  # each col holds ancestry proportions for chrs 1,...,22,X
  
  stats <- data.frame(matrix(NA,nrow=niters,ncol=2))
  colnames(stats) <- c("nonP","CAnD")
  
  for(i in 1:niters){
    # sample 23 chrs, with replacement
    chr <- sample(1:23,23,replace=TRUE)
  
    dat <- ancestProp[,chr]
  
    np <- nonParam_CAnD(dat,bonfCorr=FALSE)
    ca <- CAnD(dat,bonfCorr=FALSE)
    # want to store cand_stat = -2*sum(log(pval))
    stats$nonP[i] <- overallStatistic(np)
    stats$CAnD[i] <- overallStatistic(ca)
  }
  
  return(stats)
}
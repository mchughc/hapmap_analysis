#########################
# look at trios
# assess evidence of assortative mating

# function to resample 5000 times for each chr, x chr, and autosomes together
# calculates cor between randomly sampled mate pairs
getEst <- function(ethnData,chooseSeed,withX=TRUE,allChrs=TRUE){
  set.seed(chooseSeed)
  femList <- ethnData$coriell.id[ethnData$sex=="F"]
  malList <- ethnData$coriell.id[ethnData$sex=="M"]
  
  doOne <- function(withX,allChrs){
    
    if(withX==TRUE){
      newF <- sample(femList,length(femList),replace=FALSE)
      newM <- sample(malList,length(malList),replace=FALSE)
      chrxF <- ethnData$x[match(newF,ethnData$coriell.id)]
      chrxM <- ethnData$x[match(newM,ethnData$coriell.id)]
      cx <- cor(chrxF,chrxM)
    }else{ cx <- NA }
    
    if(allChrs==TRUE){
      newF <- sample(femList,length(femList),replace=FALSE)
      newM <- sample(malList,length(malList),replace=FALSE)
      chrAllF <- ethnData$all[match(newF,ethnData$coriell.id)]
      chrAllM <- ethnData$all[match(newM,ethnData$coriell.id)]
      ca <- cor(chrAllF,chrAllM)
    }else{ ca <- NA }
    
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
  res <- replicate(5000,doOne(withX,allChrs))
  rownames(res) <- c(1:22,"all","x")
  return(res)
}



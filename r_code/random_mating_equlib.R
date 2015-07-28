
## code to find how long until equilibrium under assortative mating

# 4 values to save each iteration
# male NAm, male Euro
# female NAm, female Euro
# maybe two matrices, one per sex, with two rows, one per ancestry

EqulibTime <- function(bigB,initFemale,initMale,autosomes=FALSE){

  ancesFemale <- data.frame(matrix(NA,nrow=2,ncol=bigB))
  ancesMale <- data.frame(matrix(NA,nrow=2,ncol=bigB))
  rownames(ancesFemale) <- rownames(ancesMale) <- c("NAm","Euro")

  ancesFemale[,1] <- initFemale
  ancesMale[,1] <- initMale
  ct <- 2
  diffMale <- 1
  
  if(!autosomes){
    while(diffMale>1e-04){
      ancesMale[,ct] <- ancesFemale[,(ct-1)]
      ancesFemale[,ct] <- (ancesMale[,(ct-1)]+ancesFemale[,(ct-1)])/2
      diffMale <- abs(ancesMale["Euro",ct]-ancesMale["Euro",(ct-1)])
      ct <- ct+1
      if(ct>bigB){diffMale <- 0}
    }
  }
  
  # make an option for autosomes too
  if(autosomes){
    while(diffMale>1e-04){
      ancesMale[,ct] <- (ancesMale[,(ct-1)]+ancesFemale[,(ct-1)])/2
      ancesFemale[,ct] <- (ancesFemale[,(ct-1)]+ancesMale[,(ct-1)])/2
      diffMale <- abs(ancesMale["Euro",ct]-ancesMale["Euro",(ct-1)])
      ct <- ct+1
      if(ct>bigB){diffMale <- 0}
    }
  }

  return(rbind(ancesMale,ancesFemale))
}

EqulibTimeBoth <- function(bigB,initFemale1,initMale1,initFemale2,initMale2){
  
  ancesFemale <- data.frame(matrix(NA,nrow=4,ncol=bigB))
  ancesMale <- data.frame(matrix(NA,nrow=4,ncol=bigB))
  rownames(ancesFemale) <- rownames(ancesMale) <- c("NAm","Euro","NAm1","Euro1")
  
  ancesFemale[,1] <- c(initFemale1,initFemale2)
  ancesMale[,1] <- c(initMale1,initMale2)
  diffMale <- 1
  
  res <- EqulibTime(2,initFemale1,initMale2)
  res2 <- EqulibTime(2,initFemale2,initMale1)
  
  ancesMale[,2] <- c(res[1:2,2],res2[1:2,2])
  ancesFemale[,2] <- c(res[3:4,2],res2[3:4,2])
  ct <- 3
  
  while(diffMale>1e-04){
    res <- EqulibTime(2,ancesFemale[1:2,(ct-1)],ancesMale[3:4,(ct-1)])
    res2 <- EqulibTime(2,ancesFemale[3:4,(ct-1)],ancesMale[1:2,(ct-1)])
    
    ancesMale[,ct] <- c(res[1:2,2],res2[1:2,2])
    ancesFemale[,ct] <- c(res[3:4,2],res2[3:4,2])
    
    diffMale <- abs(ancesMale["Euro",ct]-ancesMale["Euro",(ct-1)])
    ct <- ct+1
    if(ct>bigB){diffMale <- 0}    
  }
  rownames(ancesMale) <- paste(rownames(ancesMale),".male",sep="")
  rownames(ancesFemale) <- paste(rownames(ancesFemale),".female",sep="")
  return(rbind(ancesMale,ancesFemale))
}


initFemale <- c(0.5,0.5)
initMale <- c(1,0)
bigB <- 10000
res <- EqulibTime(bigB,initFemale,initMale)

pdf("../assortative_mating/plots/timeToEqulib_NAm.pdf")
plot(1:10,res["NAm",1:10],type="l",col="blue",xlab="Generation",ylab="Proportion NAm Ancestry",lwd=1.5)
points(1:10,res["NAm1",1:10],type="l",col="green",lwd=1.5)
legend("topright",c("Male","Female"),lty=1,col=c("blue","green"),lwd=1.5)
abline(h=2/3,lwd=0.8,col="gray")
dev.off()

initFemale <- c(0.5,0.5)
initMale <- c(0,1)
bigB <- 10000
res <- EqulibTime(bigB,initFemale,initMale)

pdf("../assortative_mating/plots/timeToEqulib_Euro.pdf")
plot(1:10,res["NAm",1:10],type="l",col="blue",xlab="Generation",ylab="Proportion NAm Ancestry",lwd=1.5)
points(1:10,res["NAm1",1:10],type="l",col="green",lwd=1.5)
legend("topright",c("Male","Female"),lty=1,col=c("blue","green"),lwd=1.5)
abline(h=1/3,lwd=0.8,col="gray")
dev.off()

initFemale <- c(1,0)
initMale <- c(0,1)
bigB <- 10000
res <- EqulibTime(bigB,initFemale,initMale)

pdf("../assortative_mating/plots/timeToEqulib_mostExtreme.pdf")
plot(1:10,res["NAm",1:10],type="l",col="blue",xlab="Generation",ylab="Proportion Ancestry",lwd=2)
points(1:10,res["NAm1",1:10],type="l",col="green",lwd=2,lty=2)
legend("topright",c("Male","Female","Autosomal"),col=c("blue","green","magenta"),lwd=2,lty=c(1,2,3))
abline(h=2/3,lwd=0.8,col="gray")
abline(h=0.5,lwd=3,col="magenta",lty=3)
dev.off()

initFemale <- c(1,0)
initMale <- c(0,1)
bigB <- 1000
res <- EqulibTime(bigB,initFemale,initMale,autosomes=TRUE)

initFemale <- c(0.25,0.75)
initMale <- c(0.5,0.5)
bigB <- 1000
res <- EqulibTime(bigB,initFemale,initMale,autosomes=TRUE)

# equlib autosomal proportions are always the mean of the female/male initial proportions for each ancestry
# always reaches equlib at first admixing generation

initFemale1 <- initMale1 <- c(0,1)
initFemale2 <- initMale2 <- c(1,0)
res <- EqulibTimeBoth(bigB,initFemale1,initMale1,initFemale2,initMale2)
# first 4 are males, last 4 are females
# first two rows are euro males, 3:4 are Nam males
# 5:6 are Euro females, 7:8 are NAm females

pdf("../assortative_mating/plots/timeToEqulib_bothMandF.pdf")
plot(1:28,res["NAm.male",1:28],type="l",col="blue",xlab="Generation",ylab="Proportion Native American Ancestry",
     lwd=1.5,main="Proportion Native American Ancestry on the X Chromosome\nAssuming Random Mating between Males and Females from each population",
     cex.main=0.8,ylim=c(0,1))
points(1:28,res["NAm1.male",1:28],type="l",col="green",lwd=1.5)
points(1:28,res["NAm2.female",1:28],type="l",col="cyan",lwd=1.5)
points(1:28,res["NAm11.female",1:28],type="l",col="magenta",lwd=1.5)
legend("topright",c("Euro Male","Euro Female","NAm Male","NAm Female"),lty=1,col=c("blue","cyan","green","magenta"),lwd=1.5)
abline(h=1/2,lwd=0.8,col="gray")
dev.off()


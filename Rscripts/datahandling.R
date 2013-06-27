## @knitr DataManipulation
emp.dists <- function(x){
  x <- cbind(x, model="empirical")
  cast.formula <- code+session+scale+enc+posold+model~result_variable   
  value.name <- "resp"
  ret.list <- list(cast(filter(x, 1), cast.formula, probabilityDist, value=value.name), cast(filter(x, 2), cast.formula, probabilityDist, value=value.name)) 
  names(ret.list) <- c("normal", "forcedguess")
  ret.list <- lapply(ret.list, rename.and.melt, 1:7, value.name)
  ret.list <- lapply(ret.list, split.sessions)
  return(ret.list) 
}
filter <- function(x, group){
  groupName <- switch(group, "normal"="normal", "forcedguess"="forcedguess")
  x.ret <- x[x$group==groupName,]
  x.ret$resp <- factor(x.ret$resp, levels=if(group==1) 0:7  else  0:5)   
  return(x.ret) 
}
probabilityDist = function(x){   
  x.freq <- table(x)
  x.relfreq <- vector(length=length(x.freq)+1)
  x.relfreq[1] <- length(x)
  # Empty cell correction:
  #x.relfreq[-1] <- (x.freq+(0.5/length(x.freq))) / (length(x)+0.5)
  # No empty Cell correction
  x.relfreq[-1] <- x.freq / length(x)
  names(x.relfreq)[1] <- "ntrials"
  names(x.relfreq)[-1] <- names(x.freq)
  return(x.relfreq) 
}
split.sessions <- function(x){
  return(list(words=x[x$session==0,], images=x[x$session==1,]))
}
rename.and.melt <- function(x, melt.id, variable.name){ 
  #  names(x)[match(paste("X",0:7,sep=""),names(x), nomatch=0)] <- as.character(0:7)
  x <- melt.data.frame(x, melt.id)
  names(x)[names(x)=="result_variable"]<-variable.name   
  return(x) 
}
print.p <- function(x){
  if(x < 0.001)
    return("< .001")
  if(x < 0.01)
    return("< .01")
  if(x < 0.05)
    return("< .05")
  if(x < 0.1)
    return(paste("= .0", round(x*100), sep=""))
  if(x <1)
    return(paste("= .", round(x*100), sep=""))
  if(x == 1)
    return("= 1.00")
  return(paste("= ", as.character(round(x, digits=2), sep="")))
}
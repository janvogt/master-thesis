filter <- function(x, group){
  groupName <- switch(group, vollenkodiert="vollenkodiert", teilenkodiert="teilenkodiert")
#  sessionID <- switch(session, first=0, second=1)
#  x.ret <- x[x$group==groupName & x$session==sessionID,]
  x.ret <- x[x$group==groupName,]
  x.ret$resp <- factor(x.ret$resp, levels=if(group==1) 0:7  else  0:5)
  return(x.ret)
}

probabilityDist = function(x){
  x.freq <- table(x)
  x.relfreq <- x.freq / length(x)
  return(x.relfreq)
}

emp.dists <- function(x){
  cast.formula <- code+session+scale+enc+posold~result_variable
  value.name <- "resp"
  ret.list <- list(cast(filter(x, 1), cast.formula, probabilityDist, value=value.name), cast(filter(x, 2), cast.formula, probabilityDist, value=value.name))
  names(ret.list) <- c("fully encoded", "partly encoded")
  return(lapply(ret.list, rename.and.melt, 1:5, value.name))
}

rename.and.melt <- function(x, melt.id, variable.name){
#  names(x)[match(paste("X",0:7,sep=""),names(x), nomatch=0)] <- as.character(0:7)
  suppressWarnings(x <- melt(x, melt.id))
  names(x)[names(x)=="result_variable"]<-variable.name
  return(x)
}

rec.data <- read.csv(file="data/Rec.csv", sep=";")

dists.list <- emp.dists(rec.data)

ggplot(subset(dists.list[[2]],session==0&code%in%levels(factor(code))), aes(x=factor(resp), y=value, group=factor(paste(enc, posold)), colour=factor(enc)))+geom_freqpoly(ymax=1, ymin=0, stat="identity")+facet_grid(code+scale~session)

filter <- function(x, group){
  groupName <- switch(group, vollenkodiert="vollenkodiert", teilenkodiert="teilenkodiert")
  x.ret <- x[x$group==groupName,]
  x.ret$resp <- factor(x.ret$resp, levels=if(group==1) 0:7  else  0:5)   
  return(x.ret) 
}

probabilityDist = function(x){   
  x.freq <- table(x)
  x.relfreq <- vector(length=length(x.freq)+1)
  x.relfreq[1] <- length(x)
  x.relfreq[-1] <- (x.freq+(0.5/length(x.freq))) / (length(x)+0.5)
  names(x.relfreq)[1] <- "ntrials"
  names(x.relfreq)[-1] <- names(x.freq)
  return(x.relfreq) 
}

emp.dists <- function(x){
  x <- cbind(x, model="empirical")
  cast.formula <- code+session+scale+enc+posold+model~result_variable   
  value.name <- "resp"
  ret.list <- list(cast(filter(x, 1), cast.formula, probabilityDist, value=value.name), cast(filter(x, 2), cast.formula, probabilityDist, value=value.name)) 
  names(ret.list) <- c("fully encoded", "partly encoded")
  ret.list <- lapply(ret.list, rename.and.melt, 1:7, value.name)
  ret.list <- lapply(ret.list, split.sessions)
  return(ret.list) 
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

plot.in.chunks <- function(data, vp.per.page){
  max.len <- length(levels(factor(data$code)))
  i <- 1
  while(i+vp.per.page-1 < max.len){
    plot.data(subset(data,code%in%levels(factor(code))[i:(i+vp.per.page-1)]))
    i <- i+vp.per.page
  }
  plot.data(subset(data,code%in%levels(factor(code))[i:max.len]))
}

plot.data <- function(data){
  print(ggplot(data,aes(x=factor(variable),y=value,group=factor(paste(enc,posold)),colour=factor(enc)))+
  geom_freqpoly(stat="identity")+
  facet_grid(code+scale~posold)+
  scale_y_continuous("Probability", limits=c(0,1)))
}

vp.data <- function(data, vp=1){
  vps <- levels(factor(data$code))
  if(!is.character(vp)) vp <- vps[vp]
  return(subset(data, code==vp))
}

model.parameter <- function(data){
  pars.grp1.ses1 <- foreach.vp(data[[1]][[1]], fit.models.to.vp.data)
  pars.grp1.ses2 <- foreach.vp(data[[1]][[2]], fit.models.to.vp.data)
  pars.grp2.ses1 <- foreach.vp(data[[2]][[1]], fit.models.to.vp.data)
  pars.grp2.ses2 <- foreach.vp(data[[2]][[2]], fit.models.to.vp.data)
  return(list(list(pars.grp1.ses1, pars.grp1.ses2), list(pars.grp2.ses1, pars.grp2.ses2)))
}

fit.models.to.vp.data <- function(vp.data){
  levels.x <- levels(factor(vp.data$variable))
  levels.enc <- levels(factor(vp.data$enc))
  levels.scales <- levels(factor(vp.data$scale))
  sdt.uv <- best.nlminb(rand.par.sdt.uv,
                        list(levels.x, levels.enc, levels.scales),
                        model.objective, 
                        model=sdt.uv, 
                        data=vp.data, 
                        lower=lower.bounds.sdt.uv(levels.x, levels.enc, levels.scales),
                        upper=upper.bounds.sdt.uv(levels.x, levels.enc, levels.scales))
  sdt.ev <- best.nlminb(rand.par.sdt.ev,
                        list(levels.x, levels.enc, levels.scales),
                        model.objective, 
                        model=sdt.ev, 
                        data=vp.data, 
                        lower=lower.bounds.sdt.ev(levels.x, levels.enc, levels.scales),
                        upper=upper.bounds.sdt.ev(levels.x, levels.enc, levels.scales))
  mpt.2.htm <- best.nlminb(rand.par.mpt.2.htm,
                           list(levels.x, levels.enc, levels.scales),
                           model.objective, 
                           model=mpt.2.htm, 
                           data=vp.data, 
                           lower=lower.bounds.mpt.2.htm(levels.x, levels.enc, levels.scales),
                           upper=upper.bounds.mpt.2.htm(levels.x, levels.enc, levels.scales))
  npar.sdt.uv <- length(lower.bounds.sdt.uv(levels.x, levels.enc, levels.scales))
  npar.sdt.ev <- length(lower.bounds.sdt.ev(levels.x, levels.enc, levels.scales))
  npar.mpt.2.htm <- length(lower.bounds.mpt.2.htm(levels.x, levels.enc, levels.scales))
  ret <- c(sdt.uv$par, sdt.uv$objective, sdt.ev$par, sdt.ev$objective, mpt.2.htm$par, mpt.2.htm$objective)
  names(ret) <- c(paste(rep("sdt.uv_",npar.sdt.uv+1), c(paste("par",1:npar.sdt.uv, sep=""), "gsq"), sep=""), 
                  paste(rep("sdt.ev_",npar.sdt.ev+1), c(paste("par",1:npar.sdt.ev, sep=""), "gsq"), sep=""), 
                  paste(rep("mpt.2.htm_",npar.mpt.2.htm), c(paste("par",1:npar.mpt.2.htm, sep=""), "gsq"), sep=""))
  return(ret)
}

best.nlminb <- function(start.fun, start.fun.par, ..., n.optim=5){
  for(run in 1:n.optim){
    if(!exists("ret")){
      ret <- nlminb(do.call(start.fun, start.fun.par), ...)
    }
    else {
      n.rep <- nlminb(do.call(start.fun, start.fun.par), ...)
      if(n.rep$objective<ret$objective) ret <- n.rep
    }
  }
  return(ret)
}

#run fun for each vp in data. Fun should return a dataframe with the same columns as data.
foreach.vp <- function(data, fun, echo=TRUE){
  ret <- NULL
  for(vp in levels(factor(data$code))){
    if(echo) cat(paste("Processing ", substitute(fun), "() for VP: ", vp, ":\n\tstart...", sep=""))
    vp.ret <- data.frame(code=vp, t(fun(vp.data(data, vp))))
    if(echo) cat("\n\tfinished.\n")
    if(is.null(ret)) ret <- vp.ret
    else ret <- rbind(ret, vp.ret)
  }
  return(ret)
}

model.objective <- function(parameter, model, data, x=variable, value=value, enc=enc, scale=scale, posold=posold, ntrials=ntrials){
  p.expected <- model(parameter, 
                      eval(substitute(x), data), 
                      eval(substitute(enc), data), 
                      eval(substitute(scale), data), 
                      eval(substitute(posold), data), 
                      levels(factor(eval(substitute(x), data))), 
                      levels(factor(eval(substitute(enc), data))), 
                      levels(factor(eval(substitute(scale), data))),
                      levels(factor(eval(substitute(posold), data))))
  n.trials <- eval(substitute(ntrials), data)
  if(!is.numeric(n.trials)) n.trials <- as.numeric(levels(n.trials)[n.trials])
  n.expected <- p.expected*as.numeric(n.trials)
  n.empirical <- eval(substitute(value), data)*n.trials
  return(2 * sum( n.empirical[n.empirical!=0] * (log(n.empirical[n.empirical!=0]) - log(n.expected[n.empirical!=0]) ) ) )
}

#x= Requested state of random variable, n.enc= number of encodings, n.scale= number of scales
#parameters: [detection(left),detection(right)]*length(enc.states[enc.states!=0]), guess(right), [detect.map*(length(x.states)/2)-1, guess.map*(length(x.states)/2)-1]*n.scale
#detect.map= 1. P(most extrem rating), 2. P(second extrem rating | not most extrem rating), 3. ...
#guess.map= 1. P(most central rating), 2. P(second central rating | not most central rating), 3. ...
rand.par.mpt.2.htm <-function(x.states, enc.states, scales){
  return(runif(2*length(enc.states[enc.states!=0])+1+(length(x.states)-2)*length(scales)))
}
lower.bounds.mpt.2.htm <- function(x.states, enc.states, scales){
  return(rep(0, 2*length(enc.states[enc.states!=0])+1+(length(x.states)-2)*length(scales)))
}
upper.bounds.mpt.2.htm <- function(x.states, enc.states, scales){
  return(rep(1, 2*length(enc.states[enc.states!=0])+1+(length(x.states)-2)*length(scales)))
}
mpt.2.htm <- function(parameter, x, enc, scale, posold, x.states, enc.states, scales, left.right.names)
{
  n.states <- length(x.states)
  if(n.states%%2 > 0) stop("x.states needs to be even.")
  enc.states <- enc.states[enc.states!=0]
  det.vec <- ifelse(enc==0, 0, parameter[(match(enc, enc.states)-1)*2+match(posold, left.right.names)])
  guess.right <- parameter[length(enc.states)*2+1]
  guess.vec <- ifelse(match(x, x.states) > n.states/2, (1-det.vec)*guess.right, (1-det.vec)*(1-guess.right))
  map.pars <- parameter[(length(enc.states)*2+2):length(parameter)]
  map.det.mat <- array(1, dim=c(length(scales), rep(n.states/2, 2)))
  map.guess.mat <- array(1, dim=c(length(scales), rep(n.states/2, 2)))
  for(s in 1:length(scales)){
    for(col in 1:(n.states/2)){
      col.par.det <- map.pars[(n.states-2) * (s-1) + (n.states/2-col) + 1]
      col.par.guess <- map.pars[(n.states/2-1) * (2*s-1) + col]
      for(row in 1:(n.states/2)){
        if(col==row & col>1) map.det.mat[s, row, col] <- col.par.det
        if(col>row) map.det.mat[s, row, col] <- 1-col.par.det
        if(col==row & col<(n.states/2)) map.guess.mat[s, row, col] <- col.par.guess
        if(col<row) map.guess.mat[s, row, col] <- 1-col.par.guess
      }
    }
  }
  det.dists <- apply(map.det.mat, c(1,2), prod)[,c((n.states/2):1, 1:(n.states/2))]
  dim(det.dists) <- c(length(scales), n.states)
  guess.dists <- apply(map.guess.mat, c(1,2), prod)[,c((n.states/2):1, 1:(n.states/2))]
  dim(guess.dists) <- c(length(scales), n.states)
  confidence.vec <- ifelse(match(posold, left.right.names, nomatch=2) == 1, (-1*match(x, x.states))+length(x.states)+1, match(x, x.states))
  possible.detect <- confidence.vec > (length(x.states)/2)
  ret <- (ifelse(possible.detect, det.vec, 0) * det.dists[(match(scale,scales)-1)*n.states + match(x, x.states)]+
          guess.vec * guess.dists[(match(scale, scales)-1)*n.states + match(x, x.states)])
  return(ret)
}

#parameters: [mu, sd]*length(enc.states[enc.states!=0]), [criterion*(length(x.states)-1)]*n.scale
#lowerBounds: c(-Inf, -Inf, ..., -Inf, 0, ...)
#upperBounds: c(Inf, Inf, ..., Inf, Inf, ...)
rand.par.sdt.uv <-function(x.states, enc.states, scales){
  rand.mu <- runif(length(enc.states[enc.states!=0]), 0.5,1.5)
  rand.sd <- runif(length(enc.states[enc.states!=0]), 1,2.5)
  rand.mem.par <- rbind(rand.mu, rand.sd)[1:(length(enc.states[enc.states!=0])*2)]
  rand.crit0 <- runif(length(scales), -1,1)
  rand.crit <- runif(length(scales)*(length(x.states)-2))
  dim(rand.crit) <- c(length(x.states)-2, length(scales))
  rand.map.par <- rbind(rand.crit0, rand.crit)[1:(length(scales)*(length(x.states)-1))]
  return(c(rand.mem.par, rand.map.par))
}
lower.bounds.sdt.uv <- function(x.states, enc.states, scales){
  return(c(rep(0:1, length(enc.states[enc.states!=0])), rep(c(-Inf, rep(0, length(x.states)-2)), length(scales))))
}
upper.bounds.sdt.uv <- function(x.states, enc.states, scales){
  return(c(rep(Inf, length(enc.states[enc.states!=0])*2 + (length(x.states)-1)*length(scales))))
}
sdt.uv <- function(parameter, x, enc, scale, posold, x.states, enc.states, scales, left.right.names){
  enc.states <- enc.states[enc.states!=0]
  mu.vec <- ifelse(enc==0, 0, parameter[(match(enc, enc.states)-1)*2+1])
  sd.vec <- ifelse(enc==0, 1, parameter[(match(enc, enc.states)-1)*2+2])
  mu.vec <- ifelse(match(posold, left.right.names, nomatch=2)==1, -mu.vec, mu.vec)
  crit <- c(NA)
  for(s in 1:length(scales)){
    offset <- length(enc.states)*2+
              (s-1)*(length(x.states)-1)
    crit <- c(crit, cumsum(parameter[(offset+1):(offset+(length(x.states)-1))]))
  }
  low.bound.vec <- ifelse(match(x, x.states)==1, -Inf, crit[(match(scale, scales)-1)*(length(x.states)-1) + match(x, x.states)])
  upper.bound.vec <- ifelse(match(x, x.states)==length(x.states), Inf, crit[1+(match(scale, scales)-1)*(length(x.states)-1) + match(x, x.states)])
  ret = pnorm(upper.bound.vec, mu.vec, sqrt(1^2+sd.vec^2)) - pnorm(low.bound.vec, mu.vec, sqrt(1^2+sd.vec^2))
  return(ret)
}

#parameters: [mu]*length(enc.states[enc.states!=0]), [criterion*(length(x.states)-1)]*n.scale
#lowerBounds: c(-Inf, -Inf, ..., -Inf, 0, ...)
#upperBounds: c(Inf, Inf, ..., Inf, Inf, ...)
rand.par.sdt.ev <-function(x.states, enc.states, scales){
  rand.mem.par <- runif(length(enc.states[enc.states!=0]), 0.5,1.5)
  rand.crit0 <- runif(length(scales), -1,1)
  rand.crit <- runif(length(scales)*(length(x.states)-2))
  dim(rand.crit) <- c(length(x.states)-2, length(scales))
  rand.map.par <- rbind(rand.crit0, rand.crit)[1:(length(scales)*(length(x.states)-1))]
  return(c(rand.mem.par, rand.map.par))
}
lower.bounds.sdt.ev <- function(x.states, enc.states, scales){
  return(c(rep(0, length(enc.states[enc.states!=0])), rep(c(-Inf, rep(0, length(x.states)-2)), length(scales))))
}
upper.bounds.sdt.ev <- function(x.states, enc.states, scales){
  return(c(rep(Inf, length(enc.states[enc.states!=0]) + (length(x.states)-1)*length(scales))))
}
sdt.ev <- function(parameter, x, enc, scale, posold, x.states, enc.states, scales, left.right.names){
  enc.states <- enc.states[enc.states!=0]
  mu.vec <- ifelse(enc==0, 0, parameter[match(enc, enc.states)])
  mu.vec <- ifelse(match(posold, left.right.names, nomatch=2)==1, -mu.vec, mu.vec)
  crit <- c(NA)
  for(s in 1:length(scales)){
    offset <- length(enc.states)+
      (s-1)*(length(x.states)-1)
    crit <- c(crit, cumsum(parameter[(offset+1):(offset+(length(x.states)-1))]))
  }
  low.bound.vec <- ifelse(match(x, x.states)==1, -Inf, crit[(match(scale, scales)-1)*(length(x.states)-1) + match(x, x.states)])
  upper.bound.vec <- ifelse(match(x, x.states)==length(x.states), Inf, crit[1+(match(scale, scales)-1)*(length(x.states)-1) + match(x, x.states)])
  ret = pnorm(upper.bound.vec, mu.vec, sqrt(1^2+1^2)) - pnorm(low.bound.vec, mu.vec, sqrt(1^2+1^2))
  return(ret)
}


suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape))
rec.data <- read.csv(file="data/Rec.csv", sep=";")
dists.list <- emp.dists(rec.data)

plot.in.chunks(dists.list[[1]][[1]], 4)
plot.in.chunks(dists.list[[1]][[2]], 8)
plot.in.chunks(dists.list[[2]][[1]], 4)
plot.in.chunks(dists.list[[2]][[2]], 8)

model.param<- model.parameter(dists.list)
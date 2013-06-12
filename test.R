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
  return(ggplot(data,aes(x=factor(variable),y=value,group=factor(paste(enc,posold)),colour=factor(enc)))+
  geom_freqpoly(stat="identity")+
  facet_grid(scale~model~posold)+
  scale_y_continuous("Probability", limits=c(0,1)))
}

plot.vp <- function(data, params){
  plot <- plot.data(data)
  mpt <-  plot.mpt.2.htm(unlist(params[grep("mpt.2.htm_par\\d+", colnames(params), perl=TRUE)]),
                          levels(factor(data$variable)), 
                          levels(factor(data$enc)), 
                          levels(factor(data$scale)))
  sdt.uv.plot <-  plot.sdt.uv(unlist(params[grep("sdt.uv_par\\d+", colnames(params), perl=TRUE)]), 
                             levels(factor(data$variable)), 
                             levels(factor(data$enc)), 
                             levels(factor(data$scale)))
  sdt.ev.plot <-  plot.sdt.ev(unlist(params[grep("sdt.ev_par\\d+", colnames(params), perl=TRUE)]), 
                             levels(factor(data$variable)), 
                             levels(factor(data$enc)), 
                             levels(factor(data$scale)))
  tmp <- ggplot_gtable(ggplot_build(sdt.ev.plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  fits <- data.frame("G^2"=round(unlist(params[grep(".*gsq", colnames(params), perl=TRUE)]), digits=2), row.names=c("UVSDT", "EVSDT", "MPT: 2HTM"))
  legend <- tmp$grobs[[leg]]
  legendArranged <- arrangeGrob(tableGrob(fits, cols="Gsq"), legend, ncol=1)
  sdt <- arrangeGrob(sdt.uv.plot+theme(legend.position="none"), sdt.ev.plot+theme(legend.position="none"), legendArranged, widths=c(3/7,3/7,1/7), nrow=1, left="UVSDT")
  models <- arrangeGrob(mpt, plot, ncol=2, widths=c(1/3,2/3))
  suppressWarnings(grid.arrange(models, sdt, ncol=1, heights=c(2/3, 1/3), main=paste("VP:", data$code[1], "Session:", data$session[1])))
}

vp.data <- function(data, vp=1){
  vps <- levels(factor(data$code))
  if(!is.character(vp)) vp <- vps[vp]
  return(subset(data, code==vp))
}


#run fun for each vp in data. Fun should return a dataframe with the same columns as data.
foreach.vp <- function(data, fun, ..., echo=TRUE){
  ret <- NULL
  further.data <- list(...)
  for(vp in levels(factor(data$code))){
    if(echo) cat(paste("Processing ", substitute(fun), "() for VP: ", vp, ":\n\tstart...", sep=""))
    if(length(further.data)>0){
      further.data.vp <- lapply(further.data, vp.data, vp)
      fun.ret <- do.call(fun, c(list(vp.data(data, vp)), further.data.vp))
      if(length(fun.ret)>0) vp.ret <- data.frame(code=vp, fun.ret)
      else vp.ret <- data.frame(code=vp, TRUE)
    }
    else{
      vp.ret <- data.frame(code=vp, t(fun(vp.data(data, vp))))
    }
    if(echo) cat("\n\tfinished.\n")
    if(is.null(ret)) ret <- vp.ret
    else ret <- rbind(ret, vp.ret)
  }
  return(ret)
}

model.parameter <- function(data){
  pars.grp1.ses1 <- foreach.vp(data[[1]][[1]], fit.models.to.vp.data)
  pars.grp1.ses2 <- foreach.vp(data[[1]][[2]], fit.models.to.vp.data)
  pars.grp2.ses1 <- foreach.vp(data[[2]][[1]], fit.models.to.vp.data)
  pars.grp2.ses2 <- foreach.vp(data[[2]][[2]], fit.models.to.vp.data)
  return(list(list(pars.grp1.ses1, pars.grp1.ses2), list(pars.grp2.ses1, pars.grp2.ses2)))
}

append.predicted.data <- function(x, param){
  x[[1]][[1]] <- foreach.vp(x[[1]][[1]], add.model.data, param[[1]][[1]])
  x[[1]][[2]] <- foreach.vp(x[[1]][[2]], add.model.data, param[[1]][[2]])
  x[[2]][[1]] <- foreach.vp(x[[2]][[1]], add.model.data, param[[2]][[1]])
  x[[2]][[2]] <- foreach.vp(x[[2]][[2]], add.model.data, param[[2]][[2]])
  return(x)
}

add.model.data <- function(x, params) {
  mpt <- data.frame(x)
  mpt$value <-  mpt.2.htm(unlist(params[grep("mpt.2.htm_par\\d+", colnames(params), perl=TRUE)]), 
                          x$variable, 
                          x$enc, 
                          x$scale, 
                          x$posold, 
                          levels(factor(x$variable)), 
                          levels(factor(x$enc)), 
                          levels(factor(x$scale)),
                          levels(factor(x$posold)))
  mpt$model <- "2HTM MPT"
  sdt.uv.df <- data.frame(x)
  sdt.uv.df$value <-  sdt.uv(unlist(params[grep("sdt.uv_par\\d+", colnames(params), perl=TRUE)]), 
                             x$variable, 
                             x$enc, 
                             x$scale, 
                             x$posold, 
                             levels(factor(x$variable)), 
                             levels(factor(x$enc)), 
                             levels(factor(x$scale)),
                             levels(factor(x$posold)))
  sdt.uv.df$model <- "UVSDT"
  sdt.ev.df <- data.frame(x)
  sdt.ev.df$value <-  sdt.ev(unlist(params[grep("sdt.ev_par\\d+", colnames(params), perl=TRUE)]), 
                             x$variable, 
                             x$enc, 
                             x$scale, 
                             x$posold, 
                             levels(factor(x$variable)), 
                             levels(factor(x$enc)), 
                             levels(factor(x$scale)),
                             levels(factor(x$posold)))
  sdt.ev.df$model <- "EVSDT"
  return(rbind(x[-1], mpt[-1], sdt.uv.df[-1], sdt.ev.df[-1]))
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
  npar.sdt.uv <- n.param.sdt.uv(levels.x, levels.enc, levels.scales)
  npar.sdt.ev <- n.param.sdt.ev(levels.x, levels.enc, levels.scales)
  npar.mpt.2.htm <- n.param.mpt.2.htm(levels.x, levels.enc, levels.scales)
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
  return(rep(0, n.param.mpt.2.htm(x.states, enc.states, scales)))
}
upper.bounds.mpt.2.htm <- function(x.states, enc.states, scales){
  return(rep(1, n.param.mpt.2.htm(x.states, enc.states, scales)))
}
n.param.mpt.2.htm <- function(x.states, enc.states, scales){
  return(2*length(enc.states[enc.states!=0])+1+(length(x.states)-2)*length(scales))
}

plot.mpt.2.htm <- function(parameter, x.states, enc.states, scales, left.right.names){
  n.states <- length(x.states)
  prob.factor <- 2
  y0.guess <- n.states/2+2
  x0 <- n.states/2+0.5
  dist.param <- parameter[(length(enc.states[enc.states!=0])*2+2):length(parameter)]
  map.dists <- mapdists.mpt.2.htm(dist.param, x.states, scales)
  for(map in scales){
    guess.dist = map.dists[[2]][match(map, levels(factor(scales))),]
    det.dist = map.dists[[1]][match(map, levels(factor(scales))),]
    param.offset <- (n.states-2)*(match(map, levels(factor(scales)))-1)
    for(mem in enc.states[enc.states!=0]){
      dl <- parameter[(match(mem, levels(factor(enc.states[enc.states!=0])))-1)*2+1]
      dr <- parameter[(match(mem, levels(factor(enc.states[enc.states!=0])))-1)*2+2]
      gr <- parameter[length(enc.states[enc.states!=0])*2+1]
      y.dist.offset <- c(rep(n.states/2+0.5, n.states), rep(y0.guess+n.states/2+0.5, n.states))
      dist.lines <- data.frame(mem=mem, 
                               map=map, 
                               x=rep(c(1:(n.states/2)-0.1, (n.states/2+1):n.states+0.1),2)-0.45, 
                               y=c(rep(n.states/2+0.5, n.states), rep(y0.guess+n.states/2+0.5, n.states)),
                               xend=rep(c(1:(n.states/2)-0.1, (n.states/2+1):n.states+0.1),2)+0.45, 
                               yend=y.dist.offset+c(det.dist, guess.dist)*prob.factor, 
                               linetype=2, 
                               color=c(rep(2,n.states), rep(3,n.states)), 
                               weight=NA, 
                               name=NA)
      axis.lines <- data.frame(mem=mem, 
                               map=map, 
                               x=c(rep(0.45,4), rep(n.states+0.55,4)),
                               y=rep(c(rep(n.states/2+0.5, 2), rep(y0.guess+n.states/2+0.5,2)), 2),
                               xend=c(rep(c(0.45, n.states/2+0.35), 2), rep(c(n.states+0.55, n.states/2+0.65), 2)), 
                               yend=rep(c(n.states/2+0.5+1*prob.factor, n.states/2+0.5, y0.guess+n.states/2+0.5+1*prob.factor, y0.guess+n.states/2+0.5), 2),
                               linetype=3, 
                               color=1, 
                               weight=NA, 
                               name=NA)
      origin <- rbind(dist.lines, axis.lines, data.frame(mem=mem, map=map, xend=c(x0+0.5,x0-0.5), yend=0, x=x0, y=y0.guess, linetype=1, color=1, weight=c(1-dr, 1-dl), name=NA))
      if(!exists("ret.val")) ret.val <- origin
      else ret.val <- rbind(ret.val, origin)
      x1.offset <- (n.states/2)/2.0
      state.lines <- data.frame(mem=mem, 
                                map=map, 
                                x=c(x0+0.5,x0-0.5, rep(x0,2)), 
                                y=c(rep(0,2),rep(y0.guess,2)), 
                                xend=c(x0+x1.offset, x0-x1.offset, 0.5+x1.offset, 0.5+n.states-x1.offset), 
                                yend=c(rep(1,2),rep(y0.guess+1,2)), 
                                linetype=1,
                                color=1,
                                weight=c(dr, dl, 1-gr, gr), 
                                name=c("d(rigth)","d(left)", "1-g(right)", "g(right)"))
      ret.val <- rbind(ret.val, state.lines)
      next.weights <- c(dr, dl, 1-gr, gr)
      for(state in 1:(n.states/2-1)){
        x0.offset <- (n.states/2-state+1)/2.0
        x1.offset <- (n.states/2-state)/2.0
        next.weights <- c(next.weights[1]*(1-dist.param[param.offset+state]),
                          next.weights[2]*(1-dist.param[param.offset+state]),
                          next.weights[3]*(1-dist.param[param.offset+n.states/2-1+state]),
                          next.weights[4]*(1-dist.param[param.offset+n.states/2-1+state]),
                          next.weights[1]*(dist.param[param.offset+state]),
                          next.weights[2]*(dist.param[param.offset+state]),
                          next.weights[3]*(dist.param[param.offset+n.states/2-1+state]),
                          next.weights[4]*(dist.param[param.offset+n.states/2-1+state]))
        map.lines <- data.frame(mem=mem, 
                                map=map, 
                                x=c(x0+x0.offset, x0-x0.offset, 0.5+x0.offset, 0.5+n.states-x0.offset), 
                                y=c(rep(state,2),rep(y0.guess+state,2)), 
                                xend=c(x0+x1.offset, x0-x1.offset, 0.5+x1.offset, 0.5+n.states-x1.offset,
                                       n.states-state+1, state, n.states/2-state+1, n.states/2+state), 
                                yend=c(rep(state+1,2),rep(y0.guess+state+1,2), rep(n.states/2, 2), rep(y0.guess+n.states/2, 2)), 
                                linetype=1,
                                color=1,
                                weight=next.weights, 
                                name=NA)
        ret.val <- rbind(ret.val, map.lines)
      }
    }
  }
  return(ggplot(data=ret.val, aes(x=xend, y=yend, xend=x, yend=y))+
          scale_colour_manual(values = c("black", "red", "blue"))+
          facet_grid(map~mem)+
          geom_text(aes(x=xend, y=y, label=ifelse(is.na(name), NA, paste(name, round(weight, digits=2), sep="=")), size=2))+
          geom_segment(aes(size=weight), lineend="round", subset=.(linetype==1))+
          geom_rect(aes(xmin=xend, xmax=x, ymin=yend, ymax=y, fill=factor(color)), subset=.(linetype==2))+
          geom_segment(aes(x=x, y=y, yend=yend, xend=xend), subset=.(linetype==3))+
          theme_bw()+
           scale_x_discrete("Response", limits=levels(factor(x.states)))+
           scale_size_continuous(range = c(0.3, 4))+
          theme(legend.position = "none", panel.background=element_blank(), axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()))
}

mapdists.mpt.2.htm <- function(map.pars, x.states, scales){
  n.states <- length(x.states)
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
  return(list(detection=det.dists, guessing=guess.dists))
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
plot.sdt.uv <- function(parameter, x.states, enc.states, scales, left.right.names){ 
  df <- data.frame()
  c.df <- data.frame()
  c.vec <- matrix(NA, ncol=length(x.states)-1, nrow=length(scales)) 
  for(scale in 1:length(scales)){
    c.vec[scale,] <- cumsum(parameter[(length(enc.states[enc.states!=0])*2+(scale-1)*(length(x.states)-1)+1):(length(enc.states[enc.states!=0])*2+(scale)*(length(x.states)-1))])
    c.df <- rbind(c.df, data.frame(scale=levels(factor(scales))[scale], cutoff=c.vec[scale,]))
  }
  x.range <- c(min(c.vec[,1])-1,
               max(c.vec[,length(x.states)-1])+1)
  x.vec <- seq(x.range[1],x.range[2], length.out=500)
  if(0 %in% enc.states) {
    mu <- 0
    sd <- sqrt(1^2+1^2)
    label <- "not shown"
    y <- dnorm(x.vec, mu, sd)
    id <- 0
    df <- rbind(df,data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu, sd=sd))
  }
  for(level in levels(factor(enc.states[enc.states!=0]))){
    mu <- parameter[(match(level, enc.states)-1)*2+1]
    sd <- parameter[(match(level, enc.states)-1)*2+2]
    attr(mu, "names")<- NULL
    attr(sd, "names")<- NULL
    label <- paste(level, "time(s) shown")
    y <- dnorm(x.vec, mu, sqrt(1+sd^2))
    id <- (match(level, levels(factor(enc.states[enc.states!=0])))-1)*2+1
    df <- rbind(df, data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu, sd=sd))
    y <- dnorm(x.vec, -mu, sqrt(1+sd^2))
    id <- (match(level, levels(factor(enc.states[enc.states!=0])))-1)*2+2
    df <- rbind(df, data.frame(id=id, x=x.vec,  y=y, label=label, mu=-mu, sd=sd))
  }
  return(ggplot(df)+facet_grid(.~scale)+geom_line(aes(x=x, y=y, colour=label, group=id), stat="identity")+geom_vline(data=c.df, aes(colour=factor(scale), xintercept=cutoff))+scale_colour_discrete("factors"))
}

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
  return(c(rep(Inf, n.param.sdt.uv(x.states, enc.states, scales))))
}
n.param.sdt.uv <- function(x.states, enc.states, scales){
  return(length(enc.states[enc.states!=0])*2 + (length(x.states)-1)*length(scales))
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
plot.sdt.ev <- function(parameter, x.states, enc.states, scales, left.right.names){
  df <- data.frame()
  c.df <- data.frame()
  c.vec <- matrix(NA, ncol=length(x.states)-1, nrow=length(scales))
  for(scale in 1:length(scales)){
    c.vec[scale,] <- cumsum(parameter[(length(enc.states[enc.states!=0])+(scale-1)*(length(x.states)-1)+1):(length(enc.states[enc.states!=0])+(scale)*(length(x.states)-1))])
    c.df <- rbind(c.df, data.frame(scale=levels(factor(scales))[scale], cutoff=c.vec[scale,]))
  }
  x.range <- c(min(c.vec[,1])-1,
               max(c.vec[,length(x.states)-1])+1)
  x.vec <- seq(x.range[1],x.range[2], length.out=500)
  if(0 %in% enc.states) {
    mu <- 0
    sd <- sqrt(1^2+1^2)
    label <- "not shown"
    y <- dnorm(x.vec, mu, sd)
    id <- 0
    df <- rbind(df,data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu, sd=sd))
  }
  for(level in levels(factor(enc.states[enc.states!=0]))){
    mu <- parameter[match(level, enc.states)]
    attr(mu, "names")<- NULL
    sd <- 1
    label <- paste(level, "time(s) shown")
    y <- dnorm(x.vec, mu, sqrt(1+sd^2))
    id <- (match(level, levels(factor(enc.states[enc.states!=0])))-1)*2+1
    df <- rbind(df, data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu, sd=sd))
    y <- dnorm(x.vec, -mu, sqrt(1+sd^2))
    id <- (match(level, levels(factor(enc.states[enc.states!=0])))-1)*2+2
    df <- rbind(df, data.frame(id=id, x=x.vec,  y=y, label=label, mu=-mu, sd=sd))
  }
  return(ggplot(df)+facet_grid(.~scale)+geom_line(aes(x=x, y=y, colour=label, group=id), stat="identity")+geom_vline(data=c.df, aes(colour=factor(scale), xintercept=cutoff))+scale_colour_discrete("factors"))
}

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
  return(c(rep(Inf, n.param.sdt.ev(x.states, enc.states, scales))))
}
n.param.sdt.ev <- function(x.states, enc.states, scales){
  return(length(enc.states[enc.states!=0]) + (length(x.states)-1)*length(scales))
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
suppressMessages(library(gridExtra))
rec.data <- read.csv(file="data/Rec.csv", sep=";")
dists.list <- emp.dists(rec.data)

#plot.in.chunks(dists.list[[1]][[1]], 1)
#plot.in.chunks(dists.list[[1]][[2]], 8)
#plot.in.chunks(dists.list[[2]][[1]], 4)
#plot.in.chunks(dists.list[[2]][[2]], 8)

model.param<- model.parameter(dists.list)
dists.list <- append.predicted.data(dists.list, model.param)


foreach.vp(dists.list[[1]][[1]], plot.vp, model.param[[1]][[1]], echo=FALSE)
foreach.vp(dists.list[[1]][[2]], plot.vp, model.param[[1]][[2]], echo=FALSE)
foreach.vp(dists.list[[2]][[1]], plot.vp, model.param[[2]][[1]], echo=FALSE)
foreach.vp(dists.list[[2]][[2]], plot.vp, model.param[[2]][[2]], echo=FALSE)
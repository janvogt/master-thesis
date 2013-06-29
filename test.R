source('Rscripts/setup.R')
source('Rscripts/datahandling.R')
source('Rscripts/fitting.R')
source('Rscripts/modelspec.R')
source('Rscripts/plotting.R')

#################
# Executed Code #
#################
rec.data <- read.csv(file="data/Rec.csv", sep=";")
dists.list <- emp.dists(rec.data)
fit.data <- do.fitting(dists.list, MPT1HTM(), MPT1HTM(TRUE), MPT1HTM(FALSE, TRUE), MPT1HTM(FALSE, TRUE, TRUE), EVSDT(), EVSDT(TRUE), EVSDT(FALSE, TRUE), EVSDT(FALSE, TRUE, TRUE))
a <- multi.gen.fit.model(MPT1HTM(), subset(dists.list[[1]][[1]], code=="AG1503"))
b <- multi.gen.fit.model(MPT1HTM(TRUE), subset(dists.list[[1]][[1]], code=="AG1503"))
c <- multi.gen.fit.model(MPT1HTM(TRUE, TRUE), subset(dists.list[[1]][[1]], code=="AG1503"))
d <- multi.gen.fit.model(MPT1HTM(FALSE, TRUE), subset(dists.list[[1]][[1]], code=="AG1503"))

a <- multi.gen.fit.model(EVSDT(), subset(dists.list[[1]][[1]], code=="AG1503"))
b <- multi.gen.fit.model(EVSDT(TRUE), subset(dists.list[[1]][[1]], code=="AG1503"))
c <- multi.gen.fit.model(EVSDT(TRUE, TRUE), subset(dists.list[[1]][[1]], code=="AG1503"))
d <- multi.gen.fit.model(EVSDT(FALSE, TRUE), subset(dists.list[[1]][[1]], code=="AG1503"))

select.mpt(list(a$mptinr, b$mptinr, c$mptinr, d$mptinr))

demographics.df <- read.csv(file="data/Demographics.csv", sep=";")
min(demographics.df$age)
mean(demographics.df$age)
table(demographics.df$occup)
table(demographics.df$occupdesc[demographics.df$occup=="student"])

plot.vp(subset(dists.list[[1]][[1]], code=="AG1503"), model.param[[1]][[1]][1, ])

########################
# Deprecated Functions #
########################
plot.UVSDT <- function(parameter, show.zero=TRUE){ 
  parameter <- parameter[grep(".*UVSDT_", names(parameter))]
  names <- names(parameter)
  #Get scales and states
  crit.scale.levels <- get.levels.from.strings(names, ".*?_c\\.(\\d+)\\.(.+)")
  scales <- crit.scale.levels[[2]]
  crit <- crit.scale.levels[[1]]
  #Get encoding strengths and left/right names
  enc.levels <- get.levels.from.strings(names, ".*?_mu\\.(.+)")
  enc.states <- enc.levels[[1]]
  left.right.names <- c("left", "right")
  #only param names
  r.res <- regexec(".*?_(.*)", names)
  names(parameter) <- mapply(function(x, y) substr(y, x[2], x[2]+attr(x, "match.length")[2]), r.res, names)
  crit.mat <- unlist(parameter[paste(rep(paste("c", crit, sep="."), each=length(scales)), scales, sep=".")])
  dim(crit.mat) <- c(length(scales), length(crit))
  crit.mat <- t(sapply(1:length(scales), function(x, y) cumsum(y[x,]), crit.mat))
  dimnames(crit.mat) <- list(scales, crit)
  mu.mat <- unlist(parameter[rep(paste("mu", enc.states, sep="."), length(left.right.names))])
  dim(mu.mat) <- c(length(enc.states), length(left.right.names))
  dimnames(mu.mat) <- list(enc.states, left.right.names)
  mu.mat[,"left"] <- -1* mu.mat[,"left"]
  sd.mat <- unlist(parameter[rep(paste("sd", enc.states, sep="."), length(left.right.names))])
  dim(sd.mat) <- c(length(enc.states), length(left.right.names))
  dimnames(sd.mat) <- list(enc.states, left.right.names)
  return(list("UVSDT"=list(mu.mat=mu.mat, sd.mat=sd.mat, crit.mat=crit.mat, show.zero=show.zero)))
}
plot.MPT2HTM <- function(parameter){
  parameter <- parameter[grep(".*MPT2HTM_", names(parameter))]
  names <- names(parameter)
  #Get scales and states
  states.scale.levels <- get.levels.from.strings(names, ".*?_mg\\.([^_]+)\\.(.+)")
  scales <- states.scale.levels[[2]]
  n.states <- (length(states.scale.levels[[1]])+1)*2  
  #Get encoding strengths and left/right names
  enc.pos.levels <- get.levels.from.strings(names, ".*?_d\\.([^_]+)\\.(.+)")
  enc.states <- enc.pos.levels[[1]]
  left.right.names <- enc.pos.levels[[2]]
  #only param names
  r.res <- regexec(".*?_(.*)", names)
  names(parameter) <-  mapply(function(x, y) substr(y, x[2], x[2]+attr(x, "match.length")[2]), r.res, names)
  dist.param <- parameter[paste(c(paste("md", 1:(n.states/2-1), sep="."), paste("mg", 1:(n.states/2-1), sep=".")), rep(scales, each=n.states-2), sep=".")]
  mapping.mat <- unlist(mapdists.mpt(dist.param, 1:n.states, scales))
  dim(mapping.mat) <- c(length(scales), n.states, 2)
  dimnames(mapping.mat) <- list(scales, 1:n.states, c("detect", "guess"))
  guessing.vec <- rep(parameter[["gr"]], length(scales))
  names(guessing.vec) <- scales
  detection.mat <- unlist(parameter[paste("d", enc.states, rep(left.right.names, each=length(enc.states)), sep=".")])
  dim(detection.mat) <- c(length(enc.states), length(left.right.names))
  dimnames(detection.mat) <- list(enc.states, left.right.names)
  return(list("MPT2HTM"=list(detection.mat=detection.mat, guessing.vec=guessing.vec, mapping.mat=mapping.mat)))
}
plot.MPT1HTM2g <- function(parameter){
  parameter <- parameter[grep(".*MPT1HTM2g_", names(parameter))]
  names <- names(parameter)
  #Get scales and states
  states.scale.levels <- get.levels.from.strings(names, ".*?_mg\\.([^_]+)\\.(.+)")
  scales <- states.scale.levels[[2]]
  n.states <- (length(states.scale.levels[[1]])+1)*2  
  #Get encoding strengths and left/right names
  enc.levels <- get.levels.from.strings(names, ".*?_d\\.(.+)")
  enc.states <- enc.levels[[1]]
  left.right.names <- c("left", "right")
  #only param names
  r.res <- regexec(".*?_(.*)", names)
  names(parameter) <-  mapply(function(x, y) substr(y, x[2], x[2]+attr(x, "match.length")[2]), r.res, names)
  dist.param <- parameter[paste(c(paste("md", 1:(n.states/2-1), sep="."), paste("mg", 1:(n.states/2-1), sep=".")), rep(scales, each=n.states-2), sep=".")]
  mapping.mat <- unlist(mapdists.mpt(dist.param, 1:n.states, scales))
  dim(mapping.mat) <- c(length(scales), n.states, 2)
  dimnames(mapping.mat) <- list(scales, 1:n.states, c("detect", "guess"))
  guessing.vec <- unlist(parameter[paste("gr", scales, sep=".")])
  names(guessing.vec) <- scales
  detection.mat <- rep(unlist(parameter[paste("d", enc.states, sep=".")]), 2)
  dim(detection.mat) <- c(length(enc.states), length(left.right.names))
  dimnames(detection.mat) <- list(enc.states, left.right.names)
  return(list("MPT1HTM2g"=list(detection.mat=detection.mat, guessing.vec=guessing.vec, mapping.mat=mapping.mat)))
}
MPT1HTM2g <- function(){
  d.vec <- function(resp, enc, scale, posold, internal=FALSE){
    enc.states <- levels(enc)[levels(enc)!=0]
    d.par <- paste("d", enc.states, sep=".")
    confidence.vec <- ifelse(level(posold)==1, -level(resp)+length(levels(resp))+1, level(resp))
    ret <- ifelse(enc==0 | (!internal & confidence.vec <= length(levels(resp))/2), NA, d.par[match(enc, enc.states)])
    return(ret)
  }
  gr.vec <- function(resp, enc, scale, posold){
    d.vec <- d.vec(resp, enc, scale, posold, internal=TRUE)
    gr.par <- paste("gr", levels(scale), sep=".")
    ret <- ifelse(level(resp)>length(levels(resp))/2, 
                  gr.par[level(scale)], 
                  comp.prob.symbol(gr.par[level(scale)]))
    ret <- ifelse(is.na(d.vec), ret, paste(comp.prob.symbol(d.vec), ret, sep="*"))
    return(ret)
  }
  dmap.vec <- function(resp, enc, scale, posold){
    dmap.par <- paste("md", rep(2:(length(levels(resp))/2)-1, length(levels(scale))), rep(levels(scale), each=length(levels(resp))/2-1), sep=".")
    dmaps <- map.par.mpt.symbol(dmap.par, length(levels(scale)))[,c(1:(length(levels(resp))/2),(length(levels(resp))/2):1)]
    ret <- dmaps[(level(resp)-1)*length(levels(scale))+level(scale)]
    return(ret)
  }
  gmap.vec <- function(resp, enc, scale, posold){
    gmap.par <- paste("mg", rep(2:(length(levels(resp))/2)-1, length(levels(scale))), rep(levels(scale), each=length(levels(resp))/2-1), sep=".")
    gmaps <- map.par.mpt.symbol(gmap.par, length(levels(scale)))[,c((length(levels(resp))/2):1,1:(length(levels(resp))/2))]
    ret <- gmaps[(level(resp)-1)*length(levels(scale))+level(scale)]
    return(ret)
  }
  model.description <- function(d.vec, gr.vec, dmap.vec, gmap.vec){
    return(ifelse(is.na(d.vec), 
                  paste(gr.vec, gmap.vec, sep="*"),
                  paste(paste(d.vec, dmap.vec, sep="*"), paste(gr.vec, gmap.vec, sep="*"), sep="+")))
  }
  fit.call <- function(data, model.filename, parameter, ...){
    return(fit.mpt(data, model.filename, ...))
  }
  return(list("MPT1HTM2g", d.vec, gr.vec, dmap.vec, gmap.vec, model.description, fit.call))
}
MPT2HTM <- function(){
  d.vec <- function(resp, enc, scale, posold, internal=FALSE){
    enc.states <- levels(enc)[levels(enc)!=0]
    d.par <- paste("d", rep(enc.states, each=length(levels(posold))), levels(posold), sep=".")
    confidence.vec <- ifelse(level(posold)==1, -level(resp)+length(levels(resp))+1, level(resp))
    ret <- ifelse(enc==0 | (!internal & confidence.vec <= length(levels(resp))/2), NA, d.par[(match(enc, enc.states)-1)*length(levels(posold))+level(posold)])
    return(ret)
  }
  gr.vec <- function(resp, enc, scale, posold){
    d.vec <- d.vec(resp, enc, scale, posold, internal=TRUE)
    gr.par <- "gr"
    ret <- ifelse(level(resp)>length(levels(resp))/2, 
                  gr.par, 
                  comp.prob.symbol(gr.par))
    ret <- ifelse(is.na(d.vec), ret, paste(comp.prob.symbol(d.vec), ret, sep="*"))
    return(ret)
  }
  dmap.vec <- function(resp, enc, scale, posold){
    dmap.par <- paste("md", rep(2:(length(levels(resp))/2)-1, length(levels(scale))), rep(levels(scale), each=length(levels(resp))/2-1), sep=".")
    dmaps <- map.par.mpt.symbol(dmap.par, length(levels(scale)))[,c(1:(length(levels(resp))/2),(length(levels(resp))/2):1)]
    ret <- dmaps[(level(resp)-1)*length(levels(scale))+level(scale)]
    return(ret)
  }
  gmap.vec <- function(resp, enc, scale, posold){
    gmap.par <- paste("mg", rep(2:(length(levels(resp))/2)-1, length(levels(scale))), rep(levels(scale), each=length(levels(resp))/2-1), sep=".")
    gmaps <- map.par.mpt.symbol(gmap.par, length(levels(scale)))[,c((length(levels(resp))/2):1,1:(length(levels(resp))/2))]
    ret <- gmaps[(level(resp)-1)*length(levels(scale))+level(scale)]
    return(ret)
  }
  model.description <- function(d.vec, gr.vec, dmap.vec, gmap.vec){
    return(ifelse(is.na(d.vec), 
                  paste(gr.vec, gmap.vec, sep="*"),
                  paste(paste(d.vec, dmap.vec, sep="*"), paste(gr.vec, gmap.vec, sep="*"), sep="+")))
  }
  fit.call <- function(data, model.filename, parameter, ...){
    return(fit.mpt(data, model.filename, ...))
  }
  return(list("MPT2HTM", d.vec, gr.vec, dmap.vec, gmap.vec, model.description, fit.call))
}
UVSDT <- function(restrictions=NULL){
  mu.vec <- function(resp, enc, scale, posold){
    enc.states <- levels(enc)[levels(enc)!=0]
    mu.par <- paste("mu", enc.states, sep=".")
    ret <- ifelse(enc==0, "0", mu.par[match(enc, enc.states)])
    ret <- ifelse(level(posold)==1, paste("-", ret, sep=""), ret)
    return(ret)
  }
  sd.vec <- function(resp, enc, scale, posold){
    enc.states <- levels(enc)[levels(enc)!=0]
    sd.par <- paste("sd", enc.states, sep=".")
    ret <- ifelse(enc==0, "1", sd.par[match(enc, enc.states)])
    return(ret)
  }
  low.bound.vec <- function(resp, enc, scale, posold){
    n.states <- length(levels(resp))
    crit.par <- sdt.crit.par(n.states, scale)
    ret <- ifelse(level(resp)==1, "-Inf", crit.par[(level(scale)-1)*(n.states) + level(resp)])
    return(ret)
  }
  upper.bound.vec <- function(resp, enc, scale, posold){
    n.states <- length(levels(resp))
    crit.par <- sdt.crit.par(n.states, scale)
    ret <- ifelse(level(resp)==n.states, "Inf", crit.par[(level(scale)-1)*(n.states) + level(resp) + 1])
    return(ret)
  }
  model.description <- function(mu.vec, sd.vec, low.bound.vec, upper.bound.vec){
    return(paste("pnorm(", upper.bound.vec, ",", mu.vec, ",sqrt(1^2+", sd.vec, "^2))-pnorm(", low.bound.vec, ",", mu.vec, ",sqrt(1^2+", sd.vec, "^2))", sep=""))
  }
  fit.call <- function(data, model.filename, parameter, ...){
    param.order <- parameter
    lower.bound <- rep(-Inf, length(param.order))
    upper.bound <- rep(Inf, length(param.order))
    lower.bound[grep("c\\.[^1].*", param.order)] <- 0
    lower.bound[grep("sd\\..*", param.order)] <- 0
    return(fit.model(data, model.filename, restrictions, use.gradient=FALSE, lower.bound=lower.bound, upper.bound=upper.bound, ...))
  }
  return(list("UVSDT", mu.vec, sd.vec, low.bound.vec, upper.bound.vec, model.description, fit.call))
}
#Accepts a model and a data.frame to fit. Returns the a data.frame with the predicted data.
gen.fit.model <- function(model.gen, data.df, ..., resp=variable, enc=enc, scale=scale, posold=posold, prob=value, model=model){
  original.data.df <- data.df
  cf <- match.call()[-2]
  cf[[1]] <- as.name("prepare.data.for.mptinr")
  data.df <- eval.parent(cf)
  model.list <- model.gen()
  par.vec.fun.list <- model.list[2:(length(model.list)-2)]
  model.func <- model.list[[length(model.list)-1]]
  fit.func <- model.list[[length(model.list)]]
  par.vec.list <- lapply(par.vec.fun.list, do.call, list(resp=data.df$resp, enc=data.df$enc, scale=data.df$scale, posold=data.df$posold))
  tmp.model <- do.call(model.func, par.vec.list)
  n.states <- length(levels(data.df$resp))
  final.model <- c()
  for(i in 1:(nrow(data.df)/n.states)){
    offset <- (i-1)*(n.states+1)
    tmp.offset <- (i-1)*n.states
    final.model[(offset+1):(offset+n.states)] <- tmp.model[(tmp.offset+1):(tmp.offset+n.states)]
    final.model[offset+n.states+1] <- ""
  }
  #print(final.model)
  mptinr.res <- fit.func(data=data.df$prob, model.filename=textConnection(final.model), parameter=check.mpt(textConnection(final.model))$parameters, ...)
  data.df[["prob"]] <- t(mptinr.res$data$predicted)
  names(data.df)[match("prob", names(data.df))] <- as.character(substitute(prob))
  ret.df <- merge(original.data.df[,-match(as.character(substitute(prob)), names(original.data.df))], data.df, by.x=c(as.character(substitute(resp)), as.character(substitute(enc)), as.character(substitute(scale)), as.character(substitute(posold))), by.y=c("resp", "enc", "scale", "posold"))
  ret.df$model <- model.list[[1]]
  return(list(data=ret.df, mptinr=mptinr.res))
}
prepare.data.for.mptinr <- function(data.df, resp=variable, enc=enc, scale=scale, posold=posold, prob=value){
  resp <- eval(substitute(resp), data.df)
  enc <- eval(substitute(enc), data.df)
  scale <- eval(substitute(scale), data.df)
  posold <- eval(substitute(posold), data.df)
  prob <- eval(substitute(prob), data.df)*data.df$ntrials
  if(!is.factor(resp)) resp <- factor(resp)
  if(!is.factor(enc)) enc <- factor(enc)
  if(!is.factor(scale)) scale <- factor(scale)
  if(!is.factor(posold)) posold <- factor(posold)
  data.df.sort <- data.frame(scale=scale, enc=enc, posold=posold, resp=resp, prob=prob)
  data.df <- data.df.sort[do.call(order, data.df.sort),]
  return(data.df)
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
  fits <- data.frame("G^2"=round(unlist(params[grep(".*gsq", colnames(params), perl=TRUE)]), digits=2), row.names=c("UVSDT", "EVSDT", "MPT: 2HTM", "MPT: 1HTM2g"))
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
    if(length(further.data)>0 & "code" %in% names(further.data)){
      further.data.vp <- lapply(further.data, vp.data, vp)
      fun.ret <- do.call(fun, c(list(vp.data(data, vp)), further.data.vp))
      if(length(fun.ret)>0) vp.ret <- data.frame(code=vp, fun.ret)
      else vp.ret <- data.frame(code=vp, TRUE)
    }
    else{
      fun.ret <- fun(vp.data(data, vp), ...)
      if(class(fun.ret)=="list"){
        fun.ret <- lapply(fun.ret, function(x, vp) return(data.frame(code=vp, x)), vp)
      }
      else{
        vp.ret <- data.frame(code=vp, t(fun(vp.data(data, vp))))
      }
    }
    if(echo) cat("\n\tfinished.\n")
    if(is.null(ret)) ret <- vp.ret
    else ret <- mapply(rbind, c(list(), ret), c(list(), vp.ret))
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
  vp.data$enc <- factor(vp.data$enc)
  mpt.1.htm.2g <- best.nlminb(rand.par.mpt.1.htm.2g,
                              list(levels.x, levels.enc, levels.scales),
                              model.objective, 
                              model=mpt.1.htm.2g, 
                              data=vp.data, 
                              lower=lower.bounds.mpt.1.htm.2g(levels.x, levels.enc, levels.scales),
                              upper=upper.bounds.mpt.1.htm.2g(levels.x, levels.enc, levels.scales))
  npar.sdt.uv <- n.param.sdt.uv(levels.x, levels.enc, levels.scales)
  npar.sdt.ev <- n.param.sdt.ev(levels.x, levels.enc, levels.scales)
  npar.mpt.2.htm <- n.param.mpt.2.htm(levels.x, levels.enc, levels.scales)
  npar.mpt.1.htm.2g <- n.param.mpt.1.htm.2g(levels.x, levels.enc, levels.scales)
  ret <- c(sdt.uv$par, sdt.uv$objective, sdt.ev$par, sdt.ev$objective, mpt.2.htm$par, mpt.2.htm$objective, mpt.1.htm.2g$par, mpt.1.htm.2g$objective)
  names(ret) <- c(paste(rep("sdt.uv_",npar.sdt.uv+1), c(paste("par",1:npar.sdt.uv, sep=""), "gsq"), sep=""), 
                  paste(rep("sdt.ev_",npar.sdt.ev+1), c(paste("par",1:npar.sdt.ev, sep=""), "gsq"), sep=""), 
                  paste(rep("mpt.2.htm_",npar.mpt.2.htm), c(paste("par",1:npar.mpt.2.htm, sep=""), "gsq"), sep=""),
                  paste(rep("mpt.1.htm.2g_",npar.mpt.1.htm.2g), c(paste("par",1:npar.mpt.1.htm.2g, sep=""), "gsq"), sep=""))
  return(ret)
}
#New (MPTinR):
gen.fit.models <- function(data, ...){
  models <- list(...)
  gen.fit.to.vp.data <- function(vp.data, models){
    fit.results <- lapply(models, gen.fit.model, vp.data[vp.data$model =="empirical",])
    predicted.data.list <- lapply(fit.results, "[[", "data")
    predicted.data <- do.call(rbind, c(predicted.data.list, list(vp.data)))
    mptinr.res.list <- lapply(fit.results, "[[", "mptinr")
    res.lines.list <- mapply(function(x, model) return(extract.mptinr.results(model()[[1]],x)), mptinr.res.list, models, SIMPLIFY=FALSE)
    res.line <- do.call(cbind, res.lines.list)
    return(list(predicted.data, res.line))
  }
  vp.data.list <- split(data, factor(data$code))
  all.vp.results <- lapply(vp.data.list, gen.fit.to.vp.data, models)
  all.data.frames <- lapply(all.vp.results, "[[", 1)
  all.fit.lines <- lapply(all.vp.results, "[[", 2)
  data.ret <- do.call(rbind, all.data.frames)
  fit.ret <- do.call(rbind, all.fit.lines)
  fit.ret <- data.frame(code=row.names(fit.ret), fit.ret)
  return(list(data=data.ret, fit.par=fit.ret))
}
extract.mptinr.results <- function(model.name, mptinr.res){
  ret <- data.frame(mptinr.res$goodness.of.fit[,1:2], mptinr.res$information.criteria[,1:2], t(mptinr.res$parameters[,1,drop=FALSE]))
  names(ret) <- paste(model.name, names(ret), sep="_")
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
  #return(2 * sum( n.empirical[n.empirical!=0] * (log(n.empirical[n.empirical!=0]) - log(n.expected[n.empirical!=0]) ) ) )
  return(-sum(eval(substitute(value), data)*log(p.expected)))
}
#x= Requested state of random variable, n.enc= number of encodings, n.scale= number of scales
#parameters: [detection]*length(enc.states[enc.states!=0]), [guess(right), detect.map*(length(x.states)/2)-1, guess.map*(length(x.states)/2)-1]*n.scale
#levels(posold) muss in richtung links-> rechts geordnet sein!
rand.par.mpt.1.htm.2g <-function(x.states, enc.states, scales){
  return(runif(n.param.mpt.1.htm.2g(x.states, enc.states, scales)))
}
lower.bounds.mpt.1.htm.2g <- function(x.states, enc.states, scales){
  return(rep(0, n.param.mpt.1.htm.2g(x.states, enc.states, scales)))
}
upper.bounds.mpt.1.htm.2g <- function(x.states, enc.states, scales){
  return(rep(1, n.param.mpt.1.htm.2g(x.states, enc.states, scales)))
}
n.param.mpt.1.htm.2g <- function(x.states, enc.states, scales){
  return(length(enc.states[enc.states!=0])+(length(x.states)-1)*length(scales))
}
mpt.1.htm.2g <- function(parameter, x, enc, scale, posold, ...){
  n.states <- length(levels(x))
  n.scales <- length(levels(scale))
  n.enc <- length(levels(enc)[levels(enc)!=0])
  if(n.states%%2 > 0) stop("Number of rating scale point has to be even.")
  det.vec <- ifelse(enc==0, 0, parameter[(match(enc, levels(enc)[levels(enc)!=0]))])
  guess.right.par <- parameter[n.enc+(1:n.scales-1)*(n.states-1)+1]
  guess.vec <- ifelse(match(x,levels(x)) > n.states/2, (1-det.vec)*guess.right.par[match(scale, levels(scale))], (1-det.vec)*(1-guess.right.par[match(scale, levels(scale))]))
  map.pars <- parameter[rep(n.enc+(1:n.scales-1)*(n.states-1), each=n.states-2)+2:(n.states-1)]
  map.dists <- mapdists.mpt(map.pars, levels(x), levels(scale))
  confidence.vec <- ifelse(match(posold, levels(posold)) == 1, (-1*match(x, levels(x)))+length(levels(x))+1, match(x, levels(x)))
  possible.detect <- confidence.vec > (length(levels(x))/2)
  ret <- (ifelse(possible.detect, det.vec, 0) * map.dists[[1]][(match(scale,levels(scale))-1)*n.states + match(x, levels(x))] +
            guess.vec * map.dists[[2]][(match(scale,levels(scale))-1)*n.states + match(x, levels(x))])
  return(ret)
}
test.mpt.1.htm.2g <- function(){
  #Test detection
  assert_equal(mpt.1.htm.2g(c(0.3, #d1
                              runif(1), #d2
                              0.45, #g1
                              runif(1), #md1
                              0.1), #mg1
                            factor(3, levels=1:4), #x
                            factor(1, levels=1:2), #enc
                            factor("low", levels=c("low")), #scale
                            factor("left", levels=c("left", "right"))),#posold
               (1-0.3)*0.45*0.1)
  assert_equal(mpt.1.htm.2g(c(runif(1), #d1
                              0.3, #d2
                              0.45, #g1
                              runif(1), #md1
                              0.1), #mg1
                            factor(3, levels=1:4), #x
                            factor(2, levels=1:2), #enc
                            factor("low", levels=c("low")), #scale
                            factor("left", levels=c("left", "right"))),#posold
               (1-0.3)*0.45*0.1)
  assert_equal(mpt.1.htm.2g(c(runif(1), #d1
                              0.3, #d2
                              0.45, #g1
                              0.2, #md1
                              0.1), #mg1
                            factor(2, levels=1:4), #x
                            factor(2, levels=1:2), #enc
                            factor("low", levels=c("low")), #scale
                            factor("left", levels=c("left", "right"))),#posold
               0.3*(1-0.2)+(1-0.3)*(1-0.45)*0.1)
  assert_equal(mpt.1.htm.2g(c(runif(1), #d1
                              0.3, #d2
                              0.45, #g1
                              0.2, #md1
                              0.1), #mg1
                            factor(4, levels=1:4), #x
                            factor(2, levels=1:2), #enc
                            factor("low", levels=c("low")), #scale
                            factor("right", levels=c("left", "right"))),#posold
               0.3*0.2+(1-0.3)*0.45*(1-0.1))
  assert_equal(mpt.1.htm.2g(c(0.7, #d1
                              runif(1), #d1
                              runif(1), #g1
                              runif(1), #md1
                              runif(1), #mg1
                              0.35, #g2
                              0.4, #md2
                              0.57), #mg2
                            factor(4, levels=1:4), #x
                            factor(1, levels=1:2), #enc
                            factor("high", levels=c("low", "high")), #scale
                            factor("right", levels=c("left", "right"))),#posold
               0.7*0.4+(1-0.7)*0.35*(1-0.57))
  assert_equal(mpt.1.htm.2g(c(0.7, #d1
                              runif(1), #d1
                              runif(1), #g1
                              runif(1), #md1
                              runif(1), #mg1
                              0.35, #g2
                              0.4, #md2
                              0.57), #mg2
                            factor(4, levels=1:4), #x
                            factor(1, levels=0:2), #enc
                            factor("high", levels=c("low", "high")), #scale
                            factor("right", levels=c("left", "right"))),#posold
               0.7*0.4+(1-0.7)*0.35*(1-0.57))
  assert_equal(mpt.1.htm.2g(c(runif(1), #d1
                              runif(1), #d1
                              runif(1), #g1
                              runif(1), #md1
                              runif(1), #mg1
                              0.35, #g2
                              0.4, #md2
                              0.57), #mg2
                            factor(4, levels=1:4), #x
                            factor(0, levels=0:2), #enc
                            factor("high", levels=c("low", "high")), #scale
                            factor("right", levels=c("left", "right"))),#posold
               0.35*(1-0.57))
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
mpt.2.htm <- function(parameter, x, enc, scale, posold, x.states, enc.states, scales, left.right.names){
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
#models: 
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
assert_equal<- function(test, expected, error=paste("Assertion Error:", substitute(test), "not equal to", substitute(expected))){
  print(test); print(expected);
  if(all(test != expected))  stop(error)
}
#plot.in.chunks(dists.list[[1]][[1]], 1)
#plot.in.chunks(dists.list[[1]][[2]], 8)
#plot.in.chunks(dists.list[[2]][[1]], 4)
#plot.in.chunks(dists.list[[2]][[2]], 8)
#fit.models.to.vp.data(dists.list[[1]][[1]][dists.list[[1]][[1]]$code=="AG1503",])
#fit.models.to.vp.data(dists.list[[2]][[1]][dists.list[[2]][[1]]$code=="RN0609",])
#model.param <- model.parameter(dists.list)
#dists.list <- append.predicted.data(dists.list, model.param)
#foreach.vp(dists.list[[1]][[1]], plot.vp, model.param[[1]][[1]], echo=FALSE)
#foreach.vp(dists.list[[1]][[2]], plot.vp, model.param[[1]][[2]], echo=FALSE)
#foreach.vp(dists.list[[2]][[1]], plot.vp, model.param[[2]][[1]], echo=FALSE)
#foreach.vp(dists.list[[2]][[2]], plot.vp, model.param[[2]][[2]], echo=FALSE)
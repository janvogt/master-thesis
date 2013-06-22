## @knitr Plotting
plot.data <- function(data){
  memory.plot <- ggplot(subset(data, enc!=0),aes(x=factor(variable),y=value,group=factor(model),colour=factor(model)))+
    geom_freqpoly(stat="identity")+
    facet_grid(enc+scale~posold)+
    scale_y_continuous("probability", limits=c(0,1))+
    scale_x_discrete("response")+
    scale_colour_discrete("model")
  if(0 %in% data$enc) {
    guessing.plot <- ggplot(subset(data, enc==0),aes(x=factor(variable),y=value,group=factor(model),colour=factor(model)))+
      geom_freqpoly(stat="identity")+
      facet_grid(posold+enc~scale)+
      scale_y_continuous("probability", limits=c(0,1))+
      scale_colour_discrete("model")
    ret.plot <- memory.plot
    tmp <- ggplot_gtable(ggplot_build(memory.plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    ret.plot <- arrangeGrob(guessing.plot+theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_text(color="white")), memory.plot+theme(legend.position="bottom"), ncol=1, heights=c(1/ (2*length(levels(factor(data$scale)))+1), 2*length(levels(factor(data$scale)))/( 2*length(levels(factor(data$scale)))+1)))
  }
  else ret.plot <- memory.plot
  return(ret.plot)
}
plot.mpt <- function(...){
  # Plots one or more MPT models
  #
  # Args:
  #   ...: one or more lists each consisting of elements:
  #           detection.mat: "numeric" matrix, dimesions: enc, left/right
  #           guessing.vec: "numeric" matrix, dimensions: scale
  #           mapping.mat: "numeric" matrix, dimensions: scale, responses, "detect"+"guessing"
  #        the argument name is taken as modelname
  #
  # Value:
  #   ggplot
  #
  arg.list <- list(...)
  .generate.df.plot.mpt <- function(model.name, detection.mat, guessing.vec, mapping.mat){
    # Plots an MPT model
    #
    # Args:
    #   detection.mat: "numeric" matrix, dimesions: enc, left/right
    #   guessing.vec: "numeric" matrix, dimensions: scale
    #   mapping.mat: "numeric" matrix, dimensions: scale, responses, "detect"+"guessing"
    #
    # Value:
    #   Dataframe with lines to draw - only for use by plot.mpt
    #
    n.states <- ncol(mapping.mat)
    prob.factor <- 2
    y0.guess <- n.states/2+2
    x0 <- n.states/2+0.5
    detection.mat <- provideDimnames(detection.mat)
    mapping.mat <- provideDimnames(mapping.mat)
    names(guessing.vec) <- dimnames(mapping.mat)[[1]]
    for(map in dimnames(mapping.mat)[[1]]){
      for(mem in dimnames(detection.mat)[[1]]){
        y.dist.offset <- c(rep(n.states/2+0.5, n.states), rep(y0.guess+n.states/2+0.5, n.states))
        mapping.weights <- c(rep(detection.mat[mem, ], each=n.states/2), rep(c(1-guessing.vec[map],guessing.vec[map]), each=n.states/2))
        dist.lines <- data.frame(mem=mem, 
                                 map=map, 
                                 x=rep(c(1:(n.states/2)-0.1, (n.states/2+1):n.states+0.1),2)-0.45, 
                                 y=c(rep(n.states/2+0.5, n.states), rep(y0.guess+n.states/2+0.5, n.states)),
                                 xend=rep(c(1:(n.states/2)-0.1, (n.states/2+1):n.states+0.1),2)+0.45, 
                                 yend=y.dist.offset+c(mapping.mat[map, , 1], mapping.mat[map, , 2])*mapping.weights*prob.factor, 
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
        origin <- rbind(dist.lines, axis.lines, data.frame(mem=mem, 
                                                           map=map, 
                                                           xend=c(x0+0.5,x0-0.5), 
                                                           yend=0, 
                                                           x=x0, 
                                                           y=y0.guess, 
                                                           linetype=1, 
                                                           color=1, 
                                                           weight=c(1-detection.mat[mem, 2], 1-detection.mat[mem, 1]), 
                                                           name=NA))
        if(!exists("ret.val")) ret.val <- origin
        else ret.val <- rbind(ret.val, origin)
        x1.offset <- (n.states/2)/2.0
        map.weights <- c(detection.mat[mem, 2], detection.mat[mem, 1], 1-guessing.vec[map], guessing.vec[map])
        state.lines <- data.frame(mem=mem, 
                                  map=map, 
                                  x=c(x0+0.5,x0-0.5, rep(x0,2)), 
                                  y=c(rep(0,2),rep(y0.guess,2)), 
                                  xend=c(x0+x1.offset, x0-x1.offset, 0.5+x1.offset, 0.5+n.states-x1.offset), 
                                  yend=c(rep(1,2),rep(y0.guess+1,2)), 
                                  linetype=1,
                                  color=1,
                                  weight=map.weights, 
                                  name=c("d(rigth)","d(left)", "1-g(right)", "g(right)"))
        ret.val <- rbind(ret.val, state.lines)
        for(state in 1:(n.states/2-1)){
          x0.offset <- (n.states/2-state+1)/2.0
          x1.offset <- (n.states/2-state)/2.0
          cur.weights <- c(map.weights[1]*sum(mapping.mat[map, (n.states/2+1):(n.states-state), 1]),
                           map.weights[2]*sum(mapping.mat[map, (state+1):(n.states/2), 1]),
                           map.weights[3]*sum(mapping.mat[map, (n.states/2+1+state):n.states, 2]),
                           map.weights[4]*sum(mapping.mat[map, 1:(n.states/2-state), 2]),
                           map.weights[1]*mapping.mat[map, n.states-state+1, 1],
                           map.weights[2]*mapping.mat[map, state, 1],
                           map.weights[3]*mapping.mat[map, n.states/2+state, 2],
                           map.weights[4]*mapping.mat[map, n.states/2-state+1, 2])
          map.lines <- data.frame(mem=mem, 
                                  map=map, 
                                  x=c(x0+x0.offset, x0-x0.offset, 0.5+x0.offset, 0.5+n.states-x0.offset), 
                                  y=c(rep(state,2),rep(y0.guess+state,2)), 
                                  xend=c(x0+x1.offset, x0-x1.offset, 0.5+x1.offset, 0.5+n.states-x1.offset,
                                         n.states-state+1, state, n.states/2-state+1, n.states/2+state), 
                                  yend=c(rep(state+1,2),rep(y0.guess+state+1,2), rep(n.states/2, 2), rep(y0.guess+n.states/2, 2)), 
                                  linetype=1,
                                  color=1,
                                  weight=cur.weights, 
                                  name=NA)
          ret.val <- rbind(ret.val, map.lines)
        }
      }
    }
    return(cbind(model=model.name, ret.val))
  }
  lines.df.list <- mapply(function(x, y) .generate.df.plot.mpt(y, x$detection.mat, x$guessing.vec, x$mapping.mat), arg.list, names(arg.list), SIMPLIFY=FALSE)
  lines.df <- do.call(rbind, lines.df.list)
  n.states <- max(sapply(arg.list, function(x) ncol(x$mapping.mat)))
  return(ggplot(data=lines.df, aes(x=xend, y=yend, xend=x, yend=y))+
           scale_colour_manual(values = c("black", "red", "blue"))+
           facet_grid(model+map~mem)+
           geom_segment(aes(size=weight), lineend="round", subset=.(linetype==1))+
           geom_segment(aes(x=x, y=y, yend=yend, xend=xend), subset=.(linetype==3))+
           geom_rect(aes(xmin=xend, xmax=x, ymin=yend, ymax=y, fill=factor(color)), subset=.(linetype==2))+
           geom_text(aes(x=xend, y=y, label=ifelse(is.na(name), NA, paste(name, round(weight, digits=2), sep="=")), size=2))+
           theme_bw()+
           scale_x_discrete("response", limits=paste("X", 1:n.states-1, sep="="))+
           scale_size_continuous(range = c(0.3, 4))+
           theme(legend.position = "none", panel.background=element_blank(), axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.line=element_blank(), panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()))
}
plot.mpts <- function(parameter){
  return(do.call(plot.mpt, c(plot.MPT2HTM(parameter), plot.MPT1HTM2g(parameter))))
}
plot.sdt <- function(...){
  # Plots one or more SDT model
  #
  # Args:
  #   ...: One named parameter for each model to be plotted. Each list must have these elements:
  #           mu.mat: "numeric" matrix, dimesions: enc, left/right
  #           sd.mat: "numeric" matrix, dimensions:  enc, left/right
  #           crit.mat: "numeric" matrix, dimensions: scale, responses
  #           show.zero: boolean, wether the not-encoded distribution should be shown or not
  #
  # Value:
  #   ggplot
  #
  arg.list <- list(...)
  .gen.df.plot.sdt <- function(model.name, mu.mat, sd.mat, x.range, show.zero=TRUE){
    mu.mat <- provideDimnames(mu.mat)
    dimnames(sd.mat) <- dimnames(mu.mat)
    df <- data.frame()
    x.vec <- seq(x.range[1],x.range[2], length.out=500)
    if(show.zero) {
      mu <- 0
      sd <- sqrt(1^2+1^2)
      label <- "not shown"
      y <- dnorm(x.vec, mu, sd)
      id <- 0
      df <- rbind(df,data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu, sd=sd))
    }
    for(enc in dimnames(mu.mat)[[1]]){
      label <- paste(enc, "time(s) shown")
      y <- dnorm(x.vec, mu.mat[enc,2], sqrt(1+sd.mat[enc,2]^2))
      id <- paste(enc, "right")
      df <- rbind(df, data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu.mat[enc,2], sd=sd.mat[enc,2]))
      y <- dnorm(x.vec, mu.mat[enc,1], sqrt(1+sd.mat[enc,1]^2))
      id <- paste(enc, "left")
      df <- rbind(df, data.frame(id=id, x=x.vec,  y=y, label=label, mu=mu.mat[enc,1], sd=sd.mat[enc,1]))
    }
    return(cbind(model=model.name, df))
  }
  .gen.crit.df.plot.sdt <- function(model.name, crit.mat){
    crit.mat <- provideDimnames(crit.mat)
    c.df <- data.frame()
    for(scale in dimnames(crit.mat)[[1]]){
      c.df <- rbind(c.df, data.frame(scale=scale, cutoff=crit.mat[scale,]))
    }
    return(cbind(model=model.name, c.df))
  }
  crit.df.list <- mapply(function(x, y) .gen.crit.df.plot.sdt(y, x$crit.mat), arg.list, names(arg.list), SIMPLIFY=FALSE)
  crit.df <- do.call(rbind, crit.df.list)
  x.range <- c(max(crit.df[["cutoff"]])*1.5, min(crit.df[["cutoff"]])*1.5)
  fun.df.list <- mapply(function(x, y, range) .gen.df.plot.sdt(y, x$mu.mat, x$sd.mat, range, x$show.zero), arg.list, names(arg.list), MoreArgs=list(x.range), SIMPLIFY=FALSE)
  fun.df <- do.call(rbind, fun.df.list)
  return(ggplot(fun.df)+facet_grid(scale~model)+geom_line(aes(x=x, y=y, colour=label, group=id), stat="identity")+geom_vline(data=crit.df, aes(xintercept=cutoff))+scale_colour_discrete("")+scale_x_continuous("familiarity")+scale_y_continuous("density"))
}
plot.sdts <- function(parameter, show.zero=TRUE){
  return(do.call(plot.sdt, c(plot.EVSDT(parameter, show.zero), plot.UVSDT(parameter, show.zero))))
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
plot.UVSDT <- function(parameter, show.zero=TRUE){ 
  arameter <- parameter[grep(".*UVSDT_", names(parameter))]
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
plot.EVSDT <- function(parameter, show.zero=TRUE){ 
  arameter <- parameter[grep(".*EVSDT_", names(parameter))]
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
  sd.mat <- matrix(1, nrow=length(enc.states), ncol=length(left.right.names))
  dimnames(sd.mat) <- list(enc.states, left.right.names)
  return(list("EVSDT"=list(mu.mat=mu.mat, sd.mat=sd.mat, crit.mat=crit.mat, show.zero=show.zero)))
}
get.levels.from.strings <- function(x, regex){
  # Extracts levels from matched subexpressions in a vector of strings.
  #
  # Args: 
  #   x: a charachter vector
  #   regex: a regular expression with subexpressions. Each of the expressions will be interpreted as a factor, which levels will be returned
  # 
  # Value:
  #   A list with one character vector for each subexpression in regex.
  match.res <- regexec(regex, x)
  matches <- regmatches(x, match.res)
  substr.list <- sapply(matches, function(x) if(length(x) > 0) return(x[2:(length(x))]))
  factors.list <- substr.list[sapply(substr.list, function(x) class(x)=="character")]
  return(lapply(1:length(factors.list[[1]]), function(x, y) levels(factor(sapply(y, "[", x))), factors.list))
}
plot.vp <- function(data, params){
  fits <- extract.fits(params)
  plot <- plot.data(data)
  if(class(plot)[1]=="gg") plot <- plot + theme(legend.position = "bottom")
  mpts <-  plot.mpts(params)
  sdts <-  plot.sdts(params, 0 %in% data$enc)+theme(legend.position = "top")
  mpt.sdt <- arrangeGrob(mpts, sdts, ncol=1, heights=c(2/3,1/3))
  models.fits <- arrangeGrob(plot, tableGrob(fits, gp=gpar(fontsize=20)), ncol=1, heights=c(5/6,1/6))
  suppressWarnings(grid.arrange(mpt.sdt, models.fits, main=paste("VP:", data$code[1], "Session:", data$session[1]+1, "Group:", ifelse(0 %in% data$enc, "'forced guessing'", "'fully encoded'")), nrow=1))
}
extract.fits <- function(params){
  models.crits.levels <- get.levels.from.strings(names(params), "([^_]*)_((?:[^cmgds]|df).*)")
  models <- models.crits.levels[[1]]
  crits <- models.crits.levels[[2]]
  crits[1:6] <- crits[c(5,4,3,6,1,2)]
  ret.mat <- unlist(params[paste(models, rep(crits, each=length(models)), sep="_")])
  dim(ret.mat) <- c(length(models), length(crits))
  dimnames(ret.mat) <- list(models, crits)
  ret.mat<- round(ret.mat, digits=3)
  ret.df <- as.data.frame(ret.mat)
  return(ret.df)
}
plot.vps <- function(data, params){
  mapply(function(x,y){
    for(vp in levels(factor(x[[1]]$code))){
      plot.vp(subset(x[[1]], code==vp), y[[1]][vp, ])
      if(vp %in% x[[2]]$code){
        plot.vp(subset(x[[2]], code==vp), y[[2]][vp, ])
      }
    }
  }, data, params)
}
mapdists.mpt <- function(map.pars, x.states, scales){
  n.states <- length(x.states)
  map.pars.vec <- unlist(map.pars)
  path.vec <- rep(1, n.states/2*2*length(scales))
  dim(map.pars.vec) <- c(n.states/2-1, 2*length(scales))
  dim(path.vec) <- c(n.states/2, 2*length(scales))
  path.vec[2:(n.states/2), ] <- sapply(data.frame(map.pars.vec), function(x) cumprod(1-x))
  path.vec[2:(n.states/2)-1, ] <- path.vec[2:(n.states/2)-1, ] * map.pars.vec
  ret.mat.det <- sapply(data.frame(path.vec[, seq(1, by=2, along.with=scales)]), function(x) return(x[c(1:length(x),length(x):1)]))
  ret.mat.guess <- sapply(data.frame(path.vec[, seq(2, by=2, along.with=scales)]), function(x) return(x[c(length(x):1,1:length(x))]))
  return(list(detection=t(ret.mat.det), guessing=t(ret.mat.guess)))
}
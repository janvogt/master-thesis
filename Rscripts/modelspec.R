## @knitr ModelSpecifications
#A model is a function returning a list of functions. 
#The first element is the model's name as single element character vector
#The second to the third to last functions are called with an equal-length factors for each: 
# function(resp, enc, scale, posold)
#and are exprected to return a same-length vector.
#The next to last function is a function(...) of all the previosly returned vectors and is expected to return a character vector of same length as its arguments representing specifying the model.
#The last function(data, model.filename, parameter, ...) is exprected to return an mptinr fit result based on it's parameters)
#Model boilerplate:
#NAME <- function(){
#   a <- function(resp, enc, scale, posold){
#     return(resp)
#   }
#   model.description <- function(a){
#     return(a)
#   }
#   fit.call <- function(data, model.filename, parameter, ...){
#     return(fit.model(data, model.filename, ...))
#   }
#   return(list("NAME", a, model.description, fit.call))
# }
MPT1HTM <- function(res.d.enc=FALSE, res.d.scale=FALSE, res.map.scale=FALSE){
  d.vec <- function(resp, enc, scale, posold, internal=FALSE){
    enc.states <- levels(enc)[levels(enc)!=0]
    d.par <- paste("d", enc.states, rep(levels(scale), each=length(enc.states)), sep=".")
    dim(d.par) <- c(length(enc.states), length(levels(scale)))
    confidence.vec <- ifelse(level(posold)==1, -level(resp)+length(levels(resp))+1, level(resp))
    ret <- ifelse(enc==0 | (!internal & confidence.vec <= length(levels(resp))/2), NA, d.par[match(enc, enc.states)+(level(scale)-1)*length(enc.states)])
    if(internal)
      return(ret)
    if(res.d.scale & length(levels(scale))>1 & res.d.enc & length(enc.states)>1)
      restrictions <- paste(d.par, collapse="=")
    else if(res.d.scale & length(levels(scale))>1)
      restrictions <- sapply(1:dim(d.par)[1], function(e,d.par) return(paste(d.par[e,], collapse="=")), d.par)
    else if(res.d.enc & length(enc.states)>1)
      restrictions <- sapply(1:dim(d.par)[2], function(e,d.par) return(paste(d.par[,e], collapse="=")), d.par)
    else
      restrictions <- NULL
    return(list(ret, restrictions))
  }
  gr.vec <- function(resp, enc, scale, posold){
    d.vec <- d.vec(resp, enc, scale, posold, internal=TRUE)
    gr.par <- paste("gr", levels(scale), sep=".")
    ret <- ifelse(level(resp)>length(levels(resp))/2, 
                  gr.par[level(scale)], 
                  comp.prob.symbol(gr.par[level(scale)]))
    ret <- ifelse(is.na(d.vec), ret, paste(comp.prob.symbol(d.vec), ret, sep="*"))
    if(res.map.scale & length(levels(scale))>1)
      restrictions <- paste(gr.par, collapse="=")
    else
      restrictions <- NULL
    return(list(ret, restrictions))
  }
  dmap.vec <- function(resp, enc, scale, posold){
    dmap.par <- paste("md", rep(2:(length(levels(resp))/2)-1, length(levels(scale))), rep(levels(scale), each=length(levels(resp))/2-1), sep=".")
    dim(dmap.par) <- c((length(levels(resp))/2)-1, length(levels(scale)))
    dmaps <- map.par.mpt.symbol(dmap.par, length(levels(scale)))[,c(1:(length(levels(resp))/2),(length(levels(resp))/2):1)]
    ret <- dmaps[(level(resp)-1)*length(levels(scale))+level(scale)]
    if(res.map.scale & length(levels(scale))>1)
      restrictions <- sapply(1:dim(dmap.par)[1], function(e,dmap.par) return(paste(dmap.par[e,], collapse="=")), dmap.par)
    else
      restrictions <- NULL
    return(list(ret, restrictions))
  }
  gmap.vec <- function(resp, enc, scale, posold){
    gmap.par <- paste("mg", rep(2:(length(levels(resp))/2)-1, length(levels(scale))), rep(levels(scale), each=length(levels(resp))/2-1), sep=".")
    dim(gmap.par) <- c((length(levels(resp))/2)-1, length(levels(scale)))
    gmaps <- map.par.mpt.symbol(gmap.par, length(levels(scale)))[,c((length(levels(resp))/2):1,1:(length(levels(resp))/2))]
    ret <- gmaps[(level(resp)-1)*length(levels(scale))+level(scale)]
    if(res.map.scale & length(levels(scale))>1)
      restrictions <- sapply(1:dim(gmap.par)[1], function(e,gmap.par) return(paste(gmap.par[e,], collapse="=")), gmap.par)
    else
      restrictions <- NULL
    return(list(ret, restrictions))
  }
  model.description <- function(d.vec, gr.vec, dmap.vec, gmap.vec){
    return(list(ifelse(is.na(d.vec[[1]]), 
                  paste(gr.vec[[1]], gmap.vec[[1]], sep="*"),
                  paste(paste(d.vec[[1]], dmap.vec[[1]], sep="*"), paste(gr.vec[[1]], gmap.vec[[1]], sep="*"), sep="+")),
                c(d.vec[[2]], gr.vec[[2]], dmap.vec[[2]], gmap.vec[[2]]))
           )
  }
  fit.call <- function(data, model.filename, restrictions, parameter, ...){
    return(fit.mpt(data, model.filename, restrictions, ...))
  }
  model.name <- ifelse(any(c(res.d.enc, res.d.scale, res.map.scale)), paste("MPT1HTM(",paste(c("res.d.enc", "res.d.scale", "res.map")[c(res.d.enc, res.d.scale, res.map.scale)], collapse=","),")", sep=""), "MPT1HTM")
  return(list(model.name, d.vec, gr.vec, dmap.vec, gmap.vec, model.description, fit.call))
}
EVSDT <- function(res.mu.enc=FALSE, res.mu.scale=FALSE, res.crit.scale=FALSE){
  mu.vec <- function(resp, enc, scale, posold){
    enc.states <- levels(enc)[levels(enc)!=0]
    mu.par <- paste("mu", enc.states, rep(levels(scale), each=length(enc.states)), sep=".")
    dim(mu.par) <- c(length(enc.states), length(levels(scale)))
    ret <- ifelse(enc==0, "0", mu.par[match(enc, enc.states)+(level(scale)-1)*length(enc.states)])
    ret <- ifelse(level(posold)==1, paste("-", ret, sep=""), ret) 
    if(res.mu.scale & length(levels(scale))>1 & res.mu.enc & length(enc.states)>1)
      restrictions <- paste(mu.par, collapse="=")
    else if(res.mu.scale & length(levels(scale))>1)
      restrictions <- sapply(1:dim(mu.par)[1], function(e,mu.par) return(paste(mu.par[e,], collapse="=")), mu.par)
    else if(res.mu.enc & length(enc.states)>1)
      restrictions <- sapply(1:dim(mu.par)[2], function(e,mu.par) return(paste(mu.par[,e], collapse="=")), mu.par)
    else
      restrictions <- NULL
    return(list(ret, restrictions))
  }
  low.bound.vec <- function(resp, enc, scale, posold){
    n.states <- length(levels(resp))
    crit.def <- sdt.crit.par(n.states, scale, res.crit.scale)
    restrictions <- crit.def[[2]]
    crit.par <- crit.def[[1]]
    ret <- ifelse(level(resp)==1, "-Inf", crit.par[(level(scale)-1)*(n.states) + level(resp)])
    return(list(ret, restrictions))
  }
  upper.bound.vec <- function(resp, enc, scale, posold){
    n.states <- length(levels(resp))
    crit.par <- sdt.crit.par(n.states, scale)[[1]]
    ret <- ifelse(level(resp)==n.states, "Inf", crit.par[(level(scale)-1)*(n.states) + level(resp) + 1])
    return(ret)
  }
  model.description <- function(mu.vec, low.bound.vec, upper.bound.vec){
    return(list(paste("pnorm((", upper.bound.vec, "-1*", mu.vec[[1]], ")/sqrt(2))-pnorm((", low.bound.vec[[1]], "-1*", mu.vec[[1]], ")/sqrt(2))", sep=""),
                c(mu.vec[[2]], low.bound.vec[[2]])))
  }
  fit.call <- function(data, model.filename, restrictions, parameter, ...){
    param.order <- parameter
    lower.bound <- rep(-Inf, length(param.order))
    upper.bound <- rep(Inf, length(param.order))
    lower.bound[grep("c\\.[^1].*", param.order)] <- 0
    return(fit.model(data, model.filename, restrictions, lower.bound=lower.bound, upper.bound=upper.bound, ...))
  }
  model.name <- ifelse(any(c(res.mu.enc, res.mu.scale, res.crit.scale)), paste("EVSDT(",paste(c("res.mu.enc", "res.mu.scale","res.crit")[c(res.mu.enc, res.mu.scale, res.crit.scale)], collapse=","),")", sep=""), "EVSDT")
  return(list(model.name, mu.vec, low.bound.vec, upper.bound.vec, model.description, fit.call))
}
map.par.mpt.symbol <- function(map.par, n.scales){
  n.states.half <- length(map.par)/n.scales+1
  maps <- array(rep("1", n.states.half*n.scales), dim=c(n.scales, n.states.half))
  for(s in 1:n.scales){
    offset <- (s-1)*(n.states.half-1)
    offset.maps <- (s-1)*(n.states.half)
    maps[s,2:n.states.half] <- cumop.symbol(comp.prob.symbol(map.par[(offset+1):(offset+n.states.half-1)]), "*", no.paren=TRUE)
    maps[s,2:n.states.half-1] <- paste(maps[s,2:n.states.half-1],
                                       map.par[(offset+1):(offset+n.states.half-1)], 
                                       sep="*")
  }
  return(maps)
}
sdt.crit.par <- function(n.states, scale, res.scale=FALSE){
  crit.par <- paste("c", rep(1:(n.states-1), length(levels(scale))), rep(levels(scale), each=n.states-1), sep=".")
  dim(crit.par) <- c(n.states-1, length(levels(scale)))
  crit <- array(NA, dim=c(n.states, length(levels(scale))))
  for(s in 1:length(levels(scale))){
    offset <- (s-1)*(n.states-1)
    crit[2:n.states, s] <- cumop.symbol(crit.par[(offset+1):(offset+n.states-1)], "+")
  }
  if(res.scale & length(levels(scale))>1)
    restrictions <- sapply(1:dim(crit.par)[1], function(e,crit.par) return(paste(crit.par[e,], collapse="=")), crit.par)
  else
    restrictions <- NULL
  return(list(crit, restrictions))
}
level <- function(factor){
  return(match(factor, levels(factor)))
}
cumop.symbol <- function(symbols, op, no.paren=FALSE){
  cumop <- c()
  for(i in 1:length(symbols)){
    tmp.cumop <- do.call(paste, c(as.list(symbols[1:i]), list(sep=op)))
    cumop[i] <- ifelse(no.paren, tmp.cumop, paste("(", tmp.cumop, ")", sep=""))
  }
  return(cumop)
}
comp.prob.symbol <- function(symbols){
  return(paste("(1-", symbols, ")", sep=""))
}
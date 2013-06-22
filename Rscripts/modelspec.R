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
EVSDT <- function(){
  mu.vec <- function(resp, enc, scale, posold){
    enc.states <- levels(enc)[levels(enc)!=0]
    mu.par <- paste("mu", enc.states, sep=".")
    ret <- ifelse(enc==0, "0", mu.par[match(enc, enc.states)])
    ret <- ifelse(level(posold)==1, paste("-", ret, sep=""), ret)
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
  model.description <- function(mu.vec, low.bound.vec, upper.bound.vec){
    return(paste("pnorm((", upper.bound.vec, "-1*", mu.vec, ")/sqrt(2))-pnorm((", low.bound.vec, "-1*", mu.vec, ")/sqrt(2))", sep=""))
  }
  fit.call <- function(data, model.filename, parameter, ...){
    param.order <- parameter
    lower.bound <- rep(-Inf, length(param.order))
    upper.bound <- rep(Inf, length(param.order))
    lower.bound[grep("c\\.[^1].*", param.order)] <- 0
    return(fit.model(data, model.filename, lower.bound=lower.bound, upper.bound=upper.bound, ...))
  }
  return(list("EVSDT", mu.vec, low.bound.vec, upper.bound.vec, model.description, fit.call))
}
UVSDT <- function(){
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
    return(fit.model(data, model.filename, use.gradient=FALSE, lower.bound=lower.bound, upper.bound=upper.bound, ...))
  }
  return(list("UVSDT", mu.vec, sd.vec, low.bound.vec, upper.bound.vec, model.description, fit.call))
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
sdt.crit.par <- function(n.states, scale){
  crit.par <- paste("c", rep(1:(n.states-1), length(levels(scale))), rep(levels(scale), each=n.states-1), sep=".")
  crit <- array(NA, dim=c(n.states, length(levels(scale))))
  for(s in 1:length(levels(scale))){
    offset <- (s-1)*(n.states-1)
    crit[2:n.states, s] <- cumop.symbol(crit.par[(offset+1):(offset+n.states-1)], "+")
  }
  return(crit)
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
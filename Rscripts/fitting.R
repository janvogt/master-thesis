## @knitr ModelFitting
do.fitting <- function(data, ...){
  ret.val.list <- lapply(data, lapply, multi.gen.fit.models, ...)
  ret.dists <- lapply(ret.val.list, lapply, "[[", "data")
  ret.dists <- mapply(function(x,y) mapply(function(x,y) rbind(x,y), x, y, SIMPLIFY = FALSE), ret.dists, data, SIMPLIFY = FALSE)
  ret.pars <- lapply(ret.val.list, lapply, "[[", "fit.par")
  ret.comp <- lapply(ret.val.list, lapply, "[[", "model.comp")
  return(list(data=ret.dists, fit.par=ret.pars, model.comp=ret.comp))
}
multi.gen.fit.models <- function(data, ...){
  models <- list(...)
  model.res.list <- lapply(models, multi.gen.fit.model, data[data$model=="empirical",])
  predicted.data.list <- lapply(model.res.list, "[[", "data")
  predicted.data.df <- do.call(rbind, predicted.data.list)
  mptinr.res.list <- lapply(model.res.list, "[[", "mptinr")
  names(mptinr.res.list) <- vapply(models, function(x) x()[[1]], "")
  fit.par.list <- mapply(function(x, name) extract.multi.mptinr.results(name, x), mptinr.res.list, names(mptinr.res.list),SIMPLIFY=FALSE, USE.NAMES=FALSE)
  fit.par <- do.call(cbind, fit.par.list)
  fit.par <- cbind(code=row.names(fit.par), fit.par)
  if(length(mptinr.res.list)>1) select.res <- select.mpt(mptinr.res.list)
  else return(list(data=predicted.data.df, fit.par=fit.par))
  return(list(data=predicted.data.df, fit.par=fit.par, model.comp=select.res))
}
multi.gen.fit.model <- function(model.gen, data.df, ..., resp=variable, enc=enc, scale=scale, posold=posold, model=model, split.by=code, prob=value, ntrials=ntrials){
  model.list <- model.gen()
  par.vec.fun.list <- model.list[2:(length(model.list)-2)]
  model.func <- model.list[[length(model.list)-1]]
  fit.func <- model.list[[length(model.list)]]
  mc <- match.call()
  mc[[1]] <- gen.model
  model.file <- eval.parent(mc[c(-4,-8)])
  mc[[1]] <- sort.data
  data.df[,as.character(substitute(value))] <- data.df[,as.character(substitute(value))]*data.df[,as.character(substitute(ntrials))]
  sorted.df.list <- by(data.df, factor(data.df[[as.character(substitute(split.by))]]), function(x, mc, model.name){
    mc[[3]] <- x[x[[model.name]]=="empirical",]
    eval.parent(mc[c(-2,-4,-8)])
  },
                       mc,
                       as.character(substitute(model)),
                       simplify=FALSE)
  data.mat <- vapply(sorted.df.list, function(x, value.name) x[[value.name]], 1:length(model.file[model.file!=""])*1.0, as.character(substitute(value)))
  mptinr.res <- fit.func(data=t(data.mat), textConnection(model.file), parameter=check.mpt(textConnection(model.file))$parameters, ...)
  ret.df <- do.call(rbind, sorted.df.list)
  predicted.data.vec <- t(mptinr.res$data$predicted$individual)
  attr(predicted.data.vec, "dim") <- NULL
  ret.df[,as.character(substitute(value))] <- predicted.data.vec/ret.df[,as.character(substitute(ntrials))]
  ret.df[,as.character(substitute(model))] <- model.list[[1]]
  row.names(ret.df) <- NULL
  return(list(data=ret.df, mptinr=mptinr.res))
}
#posold needs to be a factor with first level meaning 'left' and last level meaning 'rigth'
#returns a ready to use character vector decribing the model for mptinr
gen.model <- function(model.gen, data.df, resp=variable, enc=enc, scale=scale, posold=posold){
  mc <- match.call()
  mc[[1]] <- as.name("gen.datapoints")
  data.df <- eval.parent(mc[-2])
  model.list <- model.gen()
  par.vec.fun.list <- model.list[2:(length(model.list)-2)]
  model.func <- model.list[[length(model.list)-1]]
  par.vec.list <- lapply(par.vec.fun.list, do.call, list(resp=data.df$resp, enc=data.df$enc, scale=data.df$scale, posold=data.df$posold))
  tmp.model <- do.call(model.func, par.vec.list)
  #print(paste("Generated", model.list[[1]], ":"))
  #print(cbind(data.df, tmp.model))
  n.states <- length(levels(data.df$resp))
  final.model <- vapply(1:(nrow(data.df)/n.states), function(x, n.states, tmp.model){
    bounds <- ((x-1)*n.states+1):(x*n.states)
    return(c(tmp.model[bounds], ""))
  },
                        rep("",(n.states+1)),
                        n.states, tmp.model)
  attr(final.model, "dim") <- NULL
  return(final.model)
}
gen.datapoints <- function(data.df, resp=variable, enc=enc, scale=scale, posold=posold){
  positions <- levels(factor(eval(substitute(posold), data.df)))
  data.df <- expand.grid(resp=levels(factor(eval(substitute(resp), data.df))),
                         enc=levels(factor(eval(substitute(enc), data.df))),
                         scale=levels(factor(eval(substitute(scale), data.df))),
                         posold=positions)
  #remove are zero encodings with posold left/rigth
  data.df <- data.df[!(data.df$enc==0 & data.df$posold %in% positions[c(1,length(positions))]),]
  #remove all encoded items without position
  data.df <- data.df[!(data.df$enc!=0 & !(data.df$posold %in% positions[c(1,length(positions))])),]
  return(sort.data(data.df, resp, enc, scale, posold))
}
sort.data <- function(data.df, resp=variable, enc=enc, scale=scale, posold=posold){
  order.vec = order(data.df[[as.character(substitute(posold))]], 
                    data.df[[as.character(substitute(scale))]], 
                    data.df[[as.character(substitute(enc))]], 
                    data.df[[as.character(substitute(resp))]])
  return(data.df[order.vec,])
}
extract.multi.mptinr.results <- function(model.name, mptinr.res){
  ret <- data.frame(mptinr.res$goodness.of.fit$individual, mptinr.res$information.criteria$individual, t(mptinr.res$parameters$individual[,1,]))
  row.names(ret) <- row.names(mptinr.res$data$observed$individual)
  names(ret) <- paste(model.name, names(ret), sep="_")
  return(ret)
}
exctractModelCriteria <- function(param.mats.lists){
  ret.df <- data.frame()
  for(list.id in 1:length(param.mats.lists)){
    for(param.df.id in 1:length(param.mats.lists[[list.id]])){
      param.mat <- param.mats.lists[[list.id]][[param.df.id]]
      ret.df <- rbind(ret.df, cbind(data.frame(code=param.mat["code"], group=list.id, session=param.df.id), param.mat[,grep(".*G.Squared", names(param.mat), perl=TRUE)]))
    }
  }
  return(ret.df)
}

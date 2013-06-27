## @knitr ModelFitting
do.fitting <- function(data, ...){
  sfInit( parallel=TRUE, cpus=4)
  ret.val.list <- lapply(data, lapply, multi.gen.fit.models, ...)
  sfStop()
  ret.dists <- lapply(ret.val.list, lapply, "[[", "data")
  ret.dists <- mapply(function(x,y) mapply(function(x,y) rbind(x,y), x, y, SIMPLIFY = FALSE), ret.dists, data, SIMPLIFY = FALSE)
  ret.pars <- lapply(ret.val.list, lapply, "[[", "fit.par")
  ret.res <- lapply(ret.val.list, lapply, "[[", "model.res")
  return(list(data=ret.dists, fit.par=ret.pars, model.res=ret.res))
}
multi.gen.fit.models <- function(data, ...){
  models.lists <- list(...)
  model.res.list <- lapply(models.lists, multi.gen.fit.model, data[data$model=="empirical",], multicore=c("individual","n.optim"))
  predicted.data.list <- lapply(model.res.list, "[[", "data")
  predicted.data.df <- do.call(rbind, predicted.data.list)
  mptinr.res.list <- lapply(model.res.list, "[[", "mptinr")
  names(mptinr.res.list) <- vapply(models.lists, function(x) x[[1]], "")
  fit.par.list <- mapply(function(x, name) extract.multi.mptinr.results(name, x), mptinr.res.list, names(mptinr.res.list),SIMPLIFY=FALSE, USE.NAMES=FALSE)
  fit.par <- do.call(cbind, fit.par.list)
  fit.par <- cbind(code=row.names(fit.par), fit.par)
  return(list(data=predicted.data.df, fit.par=fit.par, model.res=mptinr.res.list))
}
multi.gen.fit.model <- function(model.list, data.df, ..., resp=variable, enc=enc, scale=scale, posold=posold, model=model, split.by=code, prob=value, ntrials=ntrials){
  par.vec.fun.list <- model.list[2:(length(model.list)-2)]
  model.func <- model.list[[length(model.list)-1]]
  fit.func <- model.list[[length(model.list)]]
  mc <- match.call()
  mc[[1]] <- gen.model
  model.def <- eval.parent(mc[c(-4,-8)])
  model.file <- model.def[[1]]
  model.restrictions <- model.def[[2]]
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
  if(length(model.restrictions)>0){
    restrictions <- as.list(model.restrictions)
    suppressWarnings(parameters <- check.mpt(textConnection(model.file), textConnection(model.restrictions))$restr.model$parameters)
  }
  else{
    restrictions <- NULL
    parameters <- check.mpt(textConnection(model.file))$parameters
  }
  mptinr.res <- fit.func(data=t(data.mat), textConnection(model.file), restrictions=restrictions, parameter=parameters, ...)
  ret.df <- do.call(rbind, sorted.df.list)
  if(is.list(mptinr.res$data$predicted)) 
    predicted.data.vec <- t(mptinr.res$data$predicted$individual) 
  else 
    predicted.data.vec <- mptinr.res$data$predicted
  attr(predicted.data.vec, "dim") <- NULL
  ret.df[,as.character(substitute(value))] <- predicted.data.vec/ret.df[,as.character(substitute(ntrials))]
  ret.df[,as.character(substitute(model))] <- model.list[[1]]
  row.names(ret.df) <- NULL
  return(list(data=ret.df, mptinr=mptinr.res))
}
#posold needs to be a factor with first level meaning 'left' and last level meaning 'rigth'
#returns a ready to use character vector decribing the model for mptinr
gen.model <- function(model.list, data.df, resp=variable, enc=enc, scale=scale, posold=posold){
  mc <- match.call()
  mc[[1]] <- as.name("gen.datapoints")
  data.df <- eval.parent(mc[-2])
  par.vec.fun.list <- model.list[2:(length(model.list)-2)]
  model.func <- model.list[[length(model.list)-1]]
  par.vec.list <- lapply(par.vec.fun.list, do.call, list(resp=data.df$resp, enc=data.df$enc, scale=data.df$scale, posold=data.df$posold))
  tmp.model <- do.call(model.func, par.vec.list)
  if(is.list(tmp.model)){
    restrictions <- tmp.model[[2]]
    tmp.model <- tmp.model[[1]]
  }
  else
    restrictions=NULL
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
  return(list(final.model, restrictions=restrictions))
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
  if(is.list(mptinr.res$parameters))
    ret <- data.frame(mptinr.res$goodness.of.fit$individual, mptinr.res$information.criteria$individual, t(mptinr.res$parameters$individual[,1,]))
  else 
    ret <- data.frame(mptinr.res$goodness.of.fit, mptinr.res$information.criteria, t(mptinr.res$parameters[,1,]))
  row.names(ret) <- row.names(mptinr.res$data$observed$individual)
  names(ret) <- paste(model.name, names(ret), sep="_")
  return(ret)
}
exctractModelCriteria <- function(param.mats.lists){
  ret.df <- data.frame()
  for(list.id in 1:length(param.mats.lists)){
    for(param.df.id in 1:length(param.mats.lists[[list.id]])){
      param.mat <- param.mats.lists[[list.id]][[param.df.id]]
      sums <- cbind(data.frame(code="(all)", group=list.id, session=param.df.id), as.list(colSums(param.mat[,grep("(.*G.Squared)|(.*df)|(.*p.value)", names(param.mat), perl=TRUE)])))
      gsqs <- unlist(sums[,grep(".*G.Squared", names(sums), perl=TRUE)])
      dfs <- unlist(sums[,grep(".*df", names(sums), perl=TRUE)])
      sums[,grep(".*p.value", names(sums), perl=TRUE)] <- as.list(1-pchisq(gsqs, dfs))
      ret.df <- rbind(ret.df, 
                      cbind(data.frame(code=param.mat["code"], group=list.id, session=param.df.id), param.mat[,grep("(.*G.Squared)|(.*df)|(.*p.value)", names(param.mat), perl=TRUE)]),
                      sums)
    }
  }
  return(ret.df)
}
gsq.test <- function(m0, m1, mptinr.res, alpha=0.05){
  if(is.list(mptinr.res[[m0]]$goodness.of.fit)){
    g.sq.diff.ind <- mptinr.res[[m1]]$goodness.of.fit$individual[["G.Squared"]]-mptinr.res[[m0]]$goodness.of.fit$individual[["G.Squared"]]
    df.diff.ind <- mptinr.res[[m1]]$goodness.of.fit$individual[["df"]]-mptinr.res[[m0]]$goodness.of.fit$individual[["df"]]
    g.sq.diff.sum <- mptinr.res[[m1]]$goodness.of.fit$sum[["G.Squared"]]-mptinr.res[[m0]]$goodness.of.fit$sum[["G.Squared"]]
    df.diff.sum <- mptinr.res[[m1]]$goodness.of.fit$sum[["df"]]-mptinr.res[[m0]]$goodness.of.fit$sum[["df"]]
    g.sq.diff.agg <- mptinr.res[[m1]]$goodness.of.fit$aggregated[["G.Squared"]]-mptinr.res[[m0]]$goodness.of.fit$aggregated[["G.Squared"]]
    df.diff.agg <- mptinr.res[[m1]]$goodness.of.fit$aggregated[["df"]]-mptinr.res[[m0]]$goodness.of.fit$aggregated[["df"]]
    g.sq.diff <- c(g.sq.diff.ind, g.sq.diff.sum, g.sq.diff.agg)
    df.diff <- c(df.diff.ind, df.diff.sum, df.diff.agg)
    codes <- c(row.names(mptinr.res[[m0]]$data$observed$individual), "sum", "aggregated")
  }
  else{
    g.sq.diff <- mptinr.res[[m1]]$goodness.of.fit[["G.Squared"]]-mptinr.res[[m0]]$goodness.of.fit[["G.Squared"]]
    df.diff <- mptinr.res[[m1]]$goodness.of.fit[["df"]]-mptinr.res[[m0]]$goodness.of.fit[["df"]]
    codes <- c(row.names(mptinr.res[[m0]]$data$observed))
  }
  p.val <- (1-pchisq(g.sq.diff,df.diff))
  sig <- p.val < alpha
  ret.df <- data.frame(code=codes, G.Squared.diff=g.sq.diff, df.diff=df.diff, p=p.val, significant=sig)
  n.sig <- if(nrow(ret.df)>1) length(which(sig[1:(length(sig)-2)]))
  ret.list <- list(ret.df, n.significant=n.sig)
  names(ret.list)[1] <- paste(m1,m0,sep="-")
  return(ret.list)
}
g.sq.diff <- function(m0, m1, mptinr.res){
  if(is.list(mptinr.res[[m0]]$goodness.of.fit)){
    g.sq.com.ind <- mptinr.res[[m1]]$goodness.of.fit$individual[["G.Squared"]]
    g.sq.base.ind <- mptinr.res[[m0]]$goodness.of.fit$individual[["G.Squared"]]
    g.sq.com.sum <- mptinr.res[[m1]]$goodness.of.fit$sum[["G.Squared"]]
    g.sq.base.sum <- mptinr.res[[m0]]$goodness.of.fit$sum[["G.Squared"]]
    g.sq.com.agg <- mptinr.res[[m1]]$goodness.of.fit$aggregated[["G.Squared"]]
    g.sq.base.agg <- mptinr.res[[m0]]$goodness.of.fit$aggregated[["G.Squared"]]
    g.sq.com <- c(g.sq.com.ind, g.sq.com.sum, g.sq.com.agg)
    g.sq.base <-  c(g.sq.base.ind, g.sq.base.sum, g.sq.base.agg)
    codes <- c(row.names(mptinr.res[[m0]]$data$observed$individual), "sum", "aggregated")
  }
  else{
    g.sq.com <- mptinr.res[[m1]]$goodness.of.fit[["G.Squared"]]
    g.sq.base <- mptinr.res[[m0]]$goodness.of.fit[["G.Squared"]]
    codes <- c(row.names(mptinr.res[[m0]]$data$observed))
  }
  g.sq.diff <- g.sq.com-g.sq.base
  who.wins <- c(m1,paste(m0,m1,sep="="),m0)[ifelse(g.sq.diff==0, 0, g.sq.diff/abs(g.sq.diff))+2]
  ret.df <- data.frame(code=codes, g.sq.base, g.sq.com, G.sq.diff=g.sq.diff, who.wins)
  names(ret.df)[c(2,3,5)] <- c(paste(m0, "G.Sq", sep="_"), paste(m1, "G.Sq", sep="_"), "superior")
  wins.df <- if(nrow(ret.df)>1){
    n.comp <- length(which(who.wins[1:(length(who.wins)-2)]==m1))
    n.base <- length(which(who.wins[1:(length(who.wins)-2)]==m0))
    n.ges <- length(who.wins)-2
    wins <- data.frame(n.base, n.comp, n.ges)
    names(wins) <- c(paste("n.wins", c(m0, m1), sep="_"), "N")
    wins
  }
  return(list(ret.df, counts=wins.df))
}
g.sq.diff.stability <- function(diff1, diff2){
  codes <- levels(factor(c(as.character(diff1[[1]]$code), as.character(diff2[[1]]$code))))
  codes <- codes[codes!="sum"&codes!="aggregated"]
  df.list <- lapply(codes, function(x){
    line1 <- diff1[[1]][diff1[[1]]$code==x,]
    line2 <- diff2[[1]][diff2[[1]]$code==x,]
    if(nrow(line1)==0)
      line1 <- data.frame(NA,NA,NA,NA,NA)
    if(nrow(line2)==0)
      line2 <- data.frame(NA,NA,NA,NA,NA)
    return(data.frame(code=x, first.ses=line1[1,4], second.ses=line2[1,4], stable=line1[1,5]==line2[1,5]))
  })
  ret.df <- do.call(rbind, df.list)
  n.first.ses <- length(which(!is.na(ret.df$first.ses)))
  n.second.ses <- length(which(!is.na(ret.df$second.ses)))
  n.both.ses <- length(which(!is.na(ret.df$stable)))
  n.stable <- length(which(ret.df$stable))
  summary <- data.frame(n.first.ses, n.second.ses, n.both.ses, n.stable)
  return(list(ret.df, summary=summary))
}

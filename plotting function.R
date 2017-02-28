library(MASS)
library(nlme)
library(formula.tools)

grps = 2
mean.list = cov.list = res.list = vector("list", grps)
mean.list[[1]] = c(10,10)
mean.list[[2]] = c(20,5)
cov.list[[1]] = matrix(c(100,0,0,25), nrow=2)
cov.list[[2]] = matrix(c(64,0,0,16), nrow=2)
res.list[[1]] = 4
res.list[[2]] = 4

mean.list 
cov.list
res.list

reff.extractor = function(model){
  model.eq = model$nlme.model
  rand.eq = model$randomFormula
  reffs = all.vars()
  predvarnames <- attr(attr(rand.eq, "terms"), "term.labels")
  
  lhs.vars(rand.eq)
}

resamp = function(model, rand.eff, res.var){
  times = seq(range(model$data$time)[1], range(model$data$time)[2], length.out=1000)
  cband = matrix(nrow=1000,ncol=2)
  resids = rnorm(nrow(rand.eff), 0, sqrt(res.var))
  for(i in 1:length(times)){
    ys = eval(getCovariateFormula(model$nlme.model)[[2]],data.frame(time=times[i],rand.eff))+resids[i]
    sort.ys = sort(ys)
    cband[i,1] = sort.ys[nrow(rand.eff)*.025+1]
    cband[i,2] = sort.ys[nrow(rand.eff)*.975]
  }
  return(cband)
}


plot.longrpart2 = function(model){
  
  #change time variable's name to "time"
  time.var = subset(rhs.vars(model$nlme.model), !(rhs.vars(model$nlme.model)%in%lhs.vars(model$fixedFormula)))
  names(model$data)[which(names(model$data)==time.var)] = "time"
  
  rand.eff = cbands = vector("list", grps)
  for(i in 1:grps){
    rand.eff[[i]] = mvrnorm(10000, mu=mean.list[[i]], Sigma=cov.list[[i]])
    colnames(rand.eff[[i]]) = lhs.vars(model$randomFormula)
   # eval(getCovariateFormula(model$nlme.model)[[2]],data.frame(time=range(model$data$time),rand.eff[[i]]))
  }
  
  
}
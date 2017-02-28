#'
#'
#' Longitudinal Recursive Partitioning Plotting Function
#'
#'
#' @param model A longrpart2 model.
#' @keywords longitudinal recursive partitioning mixed effects
#' @export
#' @examples
#' library(longRPart2)
#' # example goes here

lrp2Plot = function(model){
  # test
  #helper sub-functions
  param.extract = function(model){ #extracts parameters from model to be used in calculations
    grps = length(unique(model$leaf_node)) #number of nodes
    fix.eff = model$fixed_effects[!sapply(model$fixed_effects, is.null)] #extract fixed-effects
    res.var = model$resid.var[!sapply(model$resid.var, is.null)] #extract residual variances

    #calculate covariance matrix of random-effects from correlation matrix output by nlme
    varcov1 = model$var.corr[!sapply(model$var.corr, is.null)]
    varcov = vector("list",grps)
    for(i in 1:grps){
      varcov2 = corMatrix(varcov1[[i]])[[1]]
      varcov[[i]] = varcov2*attr(varcov2, "stdDev")*res.var[[i]]*rep(attr(varcov2, "stdDev")*res.var[[i]],each=nrow(varcov2))
      attr(varcov[[i]], "stdDev") = NULL
    }

    params = list(grps = grps, fix.eff = fix.eff, varcov = varcov, res.var = res.var)
    return(params)
  }
  mean.traj = function(model, fix.eff){ #calculates mean trajectory for plot
    times = seq(range(model$data$time)[1], range(model$data$time)[2], length.out=1000) #define time scores
    #create data frame to use in calculation
    if(model$method == "lme"){
      envir.df = data.frame(matrix(nrow=length(times),ncol=length(fix.eff)+1))
      colnames(envir.df) = c("time.scores",names(fix.eff))
      mean.est = numeric(nrow(envir.df))
      envir.df[,1] = times
      for(i in 1:nrow(envir.df)){
        envir.df[i,2:ncol(envir.df)] = fix.eff
        mean.est[i] = envir.df[i,2] + envir.df[i,3]*envir.df[i,1]
      }
      return(mean.est)
    }
    else if(model$method == "nlme"){
      envir.df = data.frame(matrix(nrow=length(times),ncol=length(fix.eff)+1))
      colnames(envir.df) = c("time",lhs.vars(model$fixedFormula))
      envir.df[,1] = times
      for(i in 1:nrow(envir.df)){
        envir.df[i,2:ncol(envir.df)] = fix.eff
      }
      mean.est = eval(getCovariateFormula(model$nlme.model)[[2]],envir.df) #evaluate model equation given times
      return(mean.est)
    }
  }
  conf.band = function(model, rand.eff, res.var){ #uses resampling to generate 95% confidence bands
    times = seq(range(model$data$time)[1], range(model$data$time)[2], length.out=1000) #define time scores
    cband = matrix(nrow=length(times),ncol=2)
    resids = rnorm(nrow(rand.eff), 0, sqrt(res.var)) #generates residuals
    #evaluates model equation a large number of times
    #calculates empirical distribution of possible scores at each time point
    #selects appropriate quantiles for 95% confidence band
    if(model$method == "lme"){
      for(i in 1:length(times)){
        ys = rand.eff[,1] + rand.eff[,2]*times[i]
        sort.ys = sort(ys)
        cband[i,1] = sort.ys[nrow(rand.eff)*.025+1]
        cband[i,2] = sort.ys[nrow(rand.eff)*.975]
      }
    }
    else if(model$method == "nlme"){
      for(i in 1:length(times)){
        ys = eval(getCovariateFormula(model$nlme.model)[[2]],data.frame(time=times[i],rand.eff))+resids[i]
        sort.ys = sort(ys)
        cband[i,1] = sort.ys[nrow(rand.eff)*.025+1]
        cband[i,2] = sort.ys[nrow(rand.eff)*.975]
      }
    }
    return(cband)
  }
  curve.struc = function(model,grps,mean.est,cbands){
    #data structure
    curve.df = data.frame(matrix(nrow=3000*grps, ncol=4))
    names(curve.df) = c("y","grp","time","col")
    #score column
    ycol=NULL
    for(i in 1:grps){
      ycol = c(ycol,mean.est[[i]],cbands[[i]][,1],cbands[[i]][,2])
    }
    curve.df$y = ycol
    #groups (to distinguish each curve)
    curve.df$grp = sort(rep(1:(3*grps),1000))
    #time variable
    times = seq(range(model$data$time)[1], range(model$data$time)[2], length.out=1000)
    curve.df$time = rep(c(rep(times,3)),grps)
    #colors
    curve.df$col = as.factor(sort(rep(1:grps, 3000)))
    return(curve.df)
  }

  #plotting function
  set.seed(1234)
  #change time variable's name to "time"
  if(model$method == "lme"){
    time.var = rhs.vars(model$fixedFormula)
    out.name = lhs.vars(model$fixedFormula)
  }
  else if(model$method == "nlme"){
    time.var = subset(rhs.vars(model$nlme.model), !(rhs.vars(model$nlme.model)%in%lhs.vars(model$fixedFormula)))
    out.name = lhs.vars(model$nlme.model)
  }
  names(model$data)[which(names(model$data)==time.var)] = "time"

  #extract group and parameter estimates
  params = param.extract(model)

  #generate mean trajectories and confidence bands
  mean.est = rand.eff = cbands = vector("list", params$grps)
  for(i in 1:params$grps){
    rand.eff[[i]] = mvrnorm(10000, mu=params$fix.eff[[i]], Sigma=params$varcov[[i]])
    colnames(rand.eff[[i]]) = lhs.vars(model$randomFormula)
    mean.est[[i]] = mean.traj(model, params$fix.eff[[i]])
    cbands[[i]] = conf.band(model,rand.eff[[i]],params$res.var[[i]])
  }
  curve.df = curve.struc(model,params$grps,mean.est,cbands)

  #plot
  p = ggplot(data=curve.df,aes(x=time, y=y, group=grp, color=col)) +
    geom_smooth(se=F) + theme_bw() +
    labs(x = time.var, y = out.name, color = "Node\n") + #
    scale_color_manual(labels = sort(unique(model$leaf_node)), values = c("blue", "red")) +
    theme(
      plot.background = element_blank()
      ,panel.grid.major = element_blank()
      ,panel.grid.minor = element_blank()
      ,panel.border = element_blank()
      ,axis.line.x = element_line(color="black")
      ,axis.line.y = element_line(color="black")
      ,legend.position="right")

  return(p)
}

#'
#'
#' Longitudinal Recursive Partitioning Plotting Function
#'
#'
#' @param model A longrpart2 model.
#' @param smooth_method Whether to use generalized additive models, smooth_method="gam",
#'                      or loess, smooth_method="loess". Defaults to loess.
#' @keywords longitudinal recursive partitioning mixed effects
#' @export
#' @examples
#' library(longRPart2)

lrp2Plot = function(model,smooth_method="loess"){

  #helper sub-functions
  param.extract = function(model) {
    grps = length(unique(model$leaf_node))
    fix.eff = model$fixed_effects[!sapply(model$fixed_effects,
                                          is.null)]
    res.var = model$resid.var[!sapply(model$resid.var, is.null)]
    varcov1 = model$var.corr[!sapply(model$var.corr, is.null)]
    varcov = vector("list", grps)
    for (i in 1:grps) {
      varcov2 = corMatrix(varcov1[[i]])[[1]]
      varcov[[i]] = varcov2 * attr(varcov2, "stdDev") *
        res.var[[i]] * rep(attr(varcov2, "stdDev") *
                             res.var[[i]], each = nrow(varcov2))
      attr(varcov[[i]], "stdDev") = NULL
    }
    params = list(grps = grps, fix.eff = fix.eff, varcov = varcov,
                  res.var = res.var)
    return(params)
  }


  mean.traj = function(model, fix.eff) {
    times = seq(range(model$data$time,na.rm=T)[1], range(model$data$time,na.rm=T)[2],
                length.out = 1000)
    if (model$method == "lme") {
      envir.df = data.frame(matrix(nrow = length(times),
                                   ncol = length(fix.eff) + 1))
      colnames(envir.df) = c("time.scores", names(fix.eff))
      mean.est = numeric(nrow(envir.df))
      envir.df[, 1] = times
      for (i in 1:nrow(envir.df)) {
        envir.df[i, 2:ncol(envir.df)] = fix.eff
        mean.est[i] = envir.df[i, 2] + envir.df[i, 3] *
          envir.df[i, 1]
      }
      return(mean.est)
    }
    else if (model$method == "nlme") {
      envir.df = data.frame(matrix(nrow = length(times),
                                   ncol = length(lhs.vars(model$fixedFormula)) + 1))
      colnames(envir.df) = c("time", lhs.vars(model$fixedFormula))
      envir.df[, 1] = times
      for (i in 1:nrow(envir.df)) {
        envir.df[i, 2:ncol(envir.df)] = fix.eff
      }

      names(envir.df)[names(envir.df) == 'time'] = time.metric

      mean.est = eval(getCovariateFormula(model$nlme.model)[[2]],
                      envir.df)
      return(mean.est)
    }
  }


  conf.band = function(model, rand.eff, res.var) {
    times = seq(range(model$data$time,na.rm=T)[1], range(model$data$time,na.rm=T)[2],
                length.out = 1000)
    cband = matrix(nrow = length(times), ncol = 2)
    resids = rnorm(nrow(rand.eff), 0, sqrt(res.var))
    if (model$method == "lme") {
      for (i in 1:length(times)) {
        ys = rand.eff[, 1] + rand.eff[, 2] * times[i]
        sort.ys = sort(ys)
        cband[i, 1] = sort.ys[nrow(rand.eff) * 0.025 +
                                1]
        cband[i, 2] = sort.ys[nrow(rand.eff) * 0.975]
      }
    }
    else if (model$method == "nlme") {
      for (i in 1:length(times)) {

        envir.for.model = data.frame(time = times[i], rand.eff)
        names(envir.for.model)[names(envir.for.model) == 'time'] = time.metric

        ys = eval(getCovariateFormula(model$nlme.model)[[2]],
                  envir.for.model) + resids[i]
        sort.ys = sort(ys)
        cband[i, 1] = sort.ys[nrow(rand.eff) * 0.025 +
                                1]
        cband[i, 2] = sort.ys[nrow(rand.eff) * 0.975]
      }
    }
    return(cband)
  }


  curve.struc = function(model, grps, mean.est, cbands) {
    curve.df = data.frame(matrix(nrow = 3000 * grps, ncol = 5))
    names(curve.df) = c("y", "grp", "time", "node", "ltype")
    ycol = NULL
    for (i in 1:grps) {
      ycol = c(ycol, mean.est[[i]], cbands[[i]][, 1], cbands[[i]][,
                                                                  2])
    }
    curve.df$y = ycol
    curve.df$grp = sort(rep(1:(3 * grps), 1000))
    curve.df$ltype = as.factor(rep(c(rep(1, 1000), rep(2,
                                                       2000)), grps))
    times = seq(range(model$data$time,na.rm=T)[1], range(model$data$time,na.rm=T)[2],
                length.out = 1000)
    curve.df$time = rep(c(rep(times, 3)), grps)
    curve.df$node = as.factor(sort(rep(unique(model$leaf_node),
                                       3000)))
    return(curve.df)
  }


  set.seed(1234)

  if (model$method == "lme") {
    time.var = rhs.vars(model$fixedFormula)
    out.name = lhs.vars(model$fixedFormula)
  }
  else if (model$method == "nlme") {
    time.var = subset(all.vars(model$nlme.model), !(all.vars(model$nlme.model) %in%
                                                      c(lhs.vars(model$nlme.model),lhs.vars(model$fixedFormula))))
    out.name = lhs.vars(model$nlme.model)


    for(k in 2:length(all.vars(model$nlme.model))){
      for(j in ncol(model$data):1){
        if(all.vars(model$nlme.model)[k] == colnames(model$data)[j]){
          time.metric = (all.vars(model$nlme.model)[k])
        }
      }
    }

    dat = model$data

    model$data$time = dat[[time.metric]]

  }


  names(model$data)[which(names(model$data) == time.var)] = "time"
  params = param.extract(model)
  mean.est = rand.eff = cbands = vector("list", params$grps)

  for (i in 1:params$grps) {

    cov.matrix = matrix(0,length(params$fix.eff[[i]]),length(params$fix.eff[[i]]))
    colnames(cov.matrix) = lhs.vars(model$fixedFormula)
    rownames(cov.matrix) = lhs.vars(model$fixedFormula)


    for(j in rownames(params$varcov[[i]])){
      for(k in colnames(params$varcov[[i]])){

        cov.matrix[j,k] = params$varcov[[i]][j,k]

      }
    }

    rand.eff[[i]] = mvrnorm(10000, mu = params$fix.eff[[i]],
                            Sigma = cov.matrix)
    colnames(rand.eff[[i]]) = lhs.vars(model$fixedFormula)
    mean.est[[i]] = mean.traj(model, params$fix.eff[[i]])
    cbands[[i]] = conf.band(model, rand.eff[[i]], params$res.var[[i]])
  }
  curve.df = curve.struc(model, params$grps, mean.est, cbands)
  model$data$node = model$leaf_node


  model$data$groupingName =   dat[[attr(terms(splitFormula(model$group, "~")[[1]]),
                                        "term.labels")]]

  model$data$responseName =   dat[[attr(terms(getResponseFormula(model$nlme.model)),
                                        "term.labels")]]

  levels(curve.df$node) = model$node2

  model$data$node <- as.factor(model$data$node)
  levels(model$data$node) <- model$node2

  p = ggplot(data = curve.df, aes(x = time, y = y, group = grp,
                                  color = "black", linetype = ltype)) + geom_smooth(se = F,method=smooth_method) +
    theme_bw() + guides(linetype = F) + facet_wrap(~node, scales = "free") +
    labs(title = paste("Node"), x = time.var, y = out.name) +
    xlim(min(curve.df$time), max(curve.df$time)) + ylim(min(c(curve.df$y,
                                                              model$data$y)),
                                                        max(c(curve.df$y, model$data$y))) +
    scale_color_manual(labels = model$node2,#sort(unique(model$leaf_node)),
                       values = c("black")) +
    theme(axis.line=element_line(),
          plot.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"), legend.position = "none") +
    geom_line(data = model$data, aes(x = time, y = responseName, group = groupingName,
                                     color = "black", linetype = NULL), alpha = 0.2) #+
  #facet_grid(~node)
  return(p)
}

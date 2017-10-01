#'
#'
#' Longitudinal Recursive Partitioning
#'
#'
#' @param method Whether to use lme() or nlme(). Use either method="lme" or method="nlme".
#'        This changes what additional arguments need to be passed.
#' @param nlme.model Necessary to specify if method="nlme"
#' @param randomFormula Random effects to include for nlme() or lme()
#' @param fixedFormula Fixed effects to include for nlme() or lme()
#' @param data Dataset
#' @param start Starting values for nlme()
#' @param group Grouping for nlme()
#' @param rPartFormula Not sure yet
#' @param weight Sample weights to be passed to rpart
#' @param R Correlation matrix to use for nlme. this is correlation=
#' @param min.dev The minimum decrease in deviance to choose a split.
#'        Note that this overrides the default cp criterion in
#'        rpart.control()
#' @param control Control function to be passed to rpart()
#' @keywords longitudinal recursive partitioning mixed effects
#' @import nlme
#' @import rpart
#' @import formula.tools
#' @import MASS
#' @import ggplot2
#' @importFrom grDevices dev.cur gray rainbow
#' @importFrom graphics legend plot points rect text
#' @importFrom stats formula lm na.omit predict quantile rnorm time nls terms
#' @export
#' @examples
#' \dontrun{
#' library(longRPart2)
#' data(ex.data.3)
#' model.0 = nlme(y~b0i+b1i*time,
#'                data=ex.data.3,
#'                fixed=b0i+b1i~1,
#'                random=b0i+b1i~1,
#'                group=~id,
#'                start=c(10,5))
#'
#'
#'lcart.mod1 <- lrp(method="nlme",
#'                  nlme.model=y~b0i+b1i*time,
#'                  fixedFormula=b0i+b1i~1,
#'                  rPartFormula = ~ z,
#'                  group= ~ id,
#'                  randomFormula=b0i+b1i~1,
#'                  data=ex.data.3,
#'                  start=c(10,5))
#'summary(lcart.mod1)
#'plot(lcart.mod1)
#'lrp2Plot(lcart.mod1)
#'}

lrp <- function(method,
                       nlme.model=NULL,
                       randomFormula,
                       fixedFormula=NULL,
                       data,
                       start,
                       group,
                       rPartFormula,
                       weight=NULL,
                       R=NULL,
                       min.dev=NULL,
                       control = rpart.control()){


  if(method=="lme"){
    lmeFormula <- fixedFormula
  }

  if(is.null(min.dev)==FALSE){

    if(method=="lme"){
      mod <- lme(lmeFormula,data=data,random=randomFormula,correlation=R,na.action=na.omit)
    }else{
      mod <- nlme(model=nlme.model,fixed=fixedFormula,data=data,
                  random=randomFormula,correlation=R,na.action=na.omit,start=start,groups=group)
    }

    control$cp <- 1 - (-2*mod$logLik - min.dev)/(-2*mod$logLik)
  }



  if(method=="lme"){
    groupingName = attr(terms(splitFormula(randomFormula,'|')[[2]]),"term.labels")
    responseName = attr(terms(getResponseFormula(lmeFormula)),"term.labels")
    groupingFactor = data[,names(data)==groupingName]
    terms = attr(terms(lmeFormula),"term.labels")
    continuous = !is.factor(data[,names(data)==terms[1]])
    ### The 3 subfunctions necessary for rpart to work.
    # The evaluation function.
    # Called once per node:
    # returns a list of two variables: a label for the node
    # and a deviance value for the node.  The deviance is
    # of length one, equal to 0 if the node is perfect/ unsplittable
    # larger equals worse
    evaluation <- function(y, wt, parms){
      model = lme(lmeFormula,data=parms[groupingFactor%in%y,],random=randomFormula,correlation=R,na.action=na.omit)
      if(continuous){
        slope = model$coefficients$fixed[2]
      }
      else{
        levels = length(levels(data[,names(data)==terms[1]]))
        y=model$coefficients$fixed[1:levels]
        x=1:levels
        slope = lm(y~x)$coefficients[2]
      }
      list(label=slope,deviance=-2*(model$logLik))
    }
  }else if(method=="nlme"){
    groupingName = attr(terms(splitFormula(group,'~')[[1]]),"term.labels")
    responseName = attr(terms(getResponseFormula(nlme.model)),"term.labels")
    groupingFactor = data[,names(data)==groupingName]

      # need to deal with categorical vs. continuous

    evaluation <- function(y, wt, parms){
      model = nlme(model=nlme.model,fixed=fixedFormula,data=parms[groupingFactor%in%y,],
                   random=randomFormula,correlation=R,na.action=na.omit,start=start,groups=group)

      # try not returning slope -- may be only for plotting
      slope=1
      list(label=slope,deviance=-2*(model$logLik))# ,pars=model$coefficients$fixed

    }

  }



  # The split function, where the work occurs.  This is used to decide
  # on the optimal split for a given covariate.
  # Called once per split candidate per node
  ### If the covariate, x, is continuous:
  # x variable is ordered
  # y is provided in the sort order of x
  # returns two vectors of length (n-1)
  #      goodness: goodness of the split, with larger numbers better. 0=no split
  #      direction: -1 = send y < cutpoint to the left
  #                  1 = send y < cutpoint to the right
  #
  ### If x is non-continuous
  # x is a set of integers (NOT the original values) defining the
  # groups for an unordered predictor.
  # Again, return goodness and direction
  #      direction: a vector of length m (number of groups), which is the applied
  #                 ordering for the unordered categories.  This is done so that
  #                 m-1 splits are performed, instead of all possible splits
  #      goodness: m-1 values, same idea as before

  ### pass in the dataset through the parms variable, with subj as y
  split <- function(y, wt, x, parms, continuous) {
    print(paste("splitting:", length(unique(x)), "values"))
    dev = vector()
    xUnique = unique(x)
    if (method == "lme") {
      rootDev = lme(lmeFormula, data = parms[groupingFactor %in%
                                               y, ], random = randomFormula, correlation = R,
                    na.action = na.omit)$logLik
    }
    else if (method == "nlme") {
      rootDev = nlme(model = nlme.model, fixed = fixedFormula,
                     data = parms[groupingFactor %in% y, ], random = randomFormula,
                     correlation = R, na.action = na.omit, start = start,
                     groups = group)$logLik
    }
    if (continuous) {
      for (i in xUnique) {
        yLeft = y[x <= i]
        yRight = y[x > i]
        if (length(yLeft) < control$minbucket || length(yRight) <
            control$minbucket) {
          dev = c(dev, 0)
        }
        else {
          if (method == "lme") {
            modelLeft = try(lme(lmeFormula, data = parms[groupingFactor %in%
                                                           yLeft, ], random = randomFormula, correlation = R,
                                na.action = na.omit), silent = TRUE)
            modelRight = try(lme(lmeFormula, data = parms[groupingFactor %in%
                                                            yRight, ], random = randomFormula, correlation = R,
                                 na.action = na.omit), silent = TRUE)
          }
          else if (method == "nlme") {
            modelLeft = try(nlme(model = nlme.model,
                                 fixed = fixedFormula, data = parms[groupingFactor %in%
                                                                      yLeft, ], random = randomFormula, correlation = R,
                                 na.action = na.omit, start = start, groups = group),
                            silent = TRUE)
            modelRight = try(nlme(model = nlme.model,
                                  fixed = fixedFormula, data = parms[groupingFactor %in%
                                                                       yRight, ], random = randomFormula, correlation = R,
                                  na.action = na.omit, start = start, groups = group),
                             silent = TRUE)
          }
          if (any(class(modelLeft) == "lme") && any(class(modelRight) ==
                                                    "lme")) {
            dev = c(dev, modelLeft$logLik + modelRight$logLik)
          }
          else {
            dev = c(dev, 0)
          }
        }
      }
      good = rep(0, length(x))
      for (i in 1:length(xUnique)) {
        good[x == xUnique[i]] = dev[i]
      }
      good = good[1:(length(good) - 1)]
      list(goodness = good + abs(rootDev) * (good != 0) *
             2, direction = rep(-1, length(good)))
    }
    else {
      order = rep(0, length(xUnique))
      response = parms[, names(parms) == responseName]
      for (i in 1:length(xUnique)) {
        order[i] = mean(response[x == xUnique[i]], na.rm = TRUE)
      }
      dir = sort(order, index.return = TRUE)$ix
      for (i in 1:(length(dir) - 1)) {
        yLeft = y[x %in% dir[1:i]]
        yRight = y[x %in% dir[(i + 1):length(dir)]]
        if (length(yLeft) < control$minbucket || length(yRight) <
            control$minbucket) {
          dev = c(dev, 0)
        }
        else {
          if (method == "lme") {
            modelLeft = try(lme(lmeFormula, data = parms[groupingFactor %in%
                                                           yLeft, ], random = randomFormula, correlation = R,
                                na.action = na.omit), silent = TRUE)
            modelRight = try(lme(lmeFormula, data = parms[groupingFactor %in%
                                                            yRight, ], random = randomFormula, correlation = R,
                                 na.action = na.omit), silent = TRUE)
          }
          else if (method == "nlme") {
            modelLeft = try(nlme(model = nlme.model,
                                 fixed = fixedFormula, data = parms[groupingFactor %in%
                                                                      yLeft, ], random = randomFormula, correlation = R,
                                 na.action = na.omit, start = start, groups = group),
                            silent = TRUE)
            modelRight = try(nlme(model = nlme.model,
                                  fixed = fixedFormula, data = parms[groupingFactor %in%
                                                                       yRight, ], random = randomFormula, correlation = R,
                                  na.action = na.omit, start = start, groups = group),
                             silent = TRUE)
          }
          if (any(class(modelLeft) == "lme") | any(class(modelLeft) ==
                                                   "nlme") && any(class(modelRight) == "lme") |
              any(class(modelRight) == "nlme")) {
            dev = c(dev, modelLeft$logLik + modelRight$logLik)
          }
          else {
            dev = c(dev, 0)
          }
        }
      }
      list(goodness = dev + abs(rootDev) * (dev != 0) *
             2, direction = dir)
    }
  }
  # The init function.  This is used, to the best of my knowledge, to initialize the process.
  # summary is used to fill print the report summary(model), and text is used to add text to
  # the plot of the tree.
  initialize <- function(y,offset,parms=0,wt){
    list(
      y=y,
      parms=parms,
      numresp=1,
      numy=1,
      summary=function(yval,dev,wt,ylevel,digits){paste("deviance (-2logLik)",format(signif(dev),3),"slope",signif(yval,2))},
      text= function(yval,dev,wt,ylevel,digits,n,use.n){
        if(!use.n){paste("m:",format(signif(yval,1)))}
        else{paste("n:",n)}
      }
    )
  }
  model <- list()
  model.rpart = rpart(paste(groupingName,c(rPartFormula)),method=list(eval=evaluation,
          split=split,init=initialize),control=control,data=data,parms=data)
  model$rpart_out <- model.rpart

  if(method=="lme"){
    model$lmeModel = lme(lmeFormula,data=data,random=randomFormula,correlation=R,na.action=na.omit)
    model$fixedFormula = lmeFormula
    model$lmeFormula = lmeFormula
  }else if(method=="nlme"){
    model$nlmeModel <- nlme(model=nlme.model,fixed=fixedFormula,data=data,
         random=randomFormula,correlation=R,na.action=na.omit,start=start,groups=group)
    model$fixedFormula <- fixedFormula
    model$nlme.model <- nlme.model
  }

  model$leaf_node <- model.rpart$where
  summary = fixed_effects = var.corr = resid.var = list()
  for(j in 1:length(table(model.rpart$where))){
    id <- names(table(model.rpart$where))[j]==model.rpart$where

    if(method=="lme"){
      model.out = lme(lmeFormula,data=data[id,],random=randomFormula,correlation=R,na.action=na.omit)
      summary[[as.numeric(names(table(model.rpart$where)))[j]]] <- summary(model.out)
      fixed_effects[[as.numeric(names(table(model.rpart$where)))[j]]] <- fixed.effects(model.out)
      var.corr[[as.numeric(names(table(model.rpart$where)))[j]]] <-model.out$modelStruct[[1]]
      resid.var[[as.numeric(names(table(model.rpart$where)))[j]]] <- model.out$sigma
    }else if(method=="nlme"){
      model.out <- nlme(model=nlme.model,fixed=fixedFormula,data=data[id,],
                              random=randomFormula,correlation=R,na.action=na.omit,start=start,groups=group)
      summary[[as.numeric(names(table(model.rpart$where)))[j]]] <- summary(model.out)
      fixed_effects[[as.numeric(names(table(model.rpart$where)))[j]]] <- fixed.effects(model.out)
      var.corr[[as.numeric(names(table(model.rpart$where)))[j]]] <- model.out$modelStruct[[1]]
      resid.var[[as.numeric(names(table(model.rpart$where)))[j]]] <- model.out$sigma
    }
  }

  model$summary <- summary
  model$fixed_effects <- fixed_effects
  model$var.corr <- var.corr
  model$resid.var <- resid.var
  model$rpart_out <- model.rpart
  model$randomFormula = randomFormula
  model$R = R

  if(method=="nlme"){
    model$group = group
    model$start = start
  }
  model$method = method
  model$data = data
  model$groupingName = groupingName
  model$rPartFormula = rPartFormula

  model$call <- match.call()
  class(model) <- "lrp"

  model
}

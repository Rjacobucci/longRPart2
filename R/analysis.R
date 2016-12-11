#
# trying to cross-validate
#
lrpCV <- function(model){
#  set.seed(2)
  # extracting the rpart-values
  dat <- model$data
  control <- model$control
  evaluation = model$functions$eval
  split = model$functions$split
  initialize = model$functions$init
  groupingName = model$groupingName
  rPartFormula = model$rPartFormula
  v <- model$control$xval
  if(v < 2)
    stop("cannot cross-validate with fewer than 2 folds.")
  n = length(unique(dat[,groupingName]))
  deviances <- data.frame(matrix(nrow=n,ncol=2))
  names(deviances) <- c(groupingName,'deviance')
  z = 1
  subj = sample(unique(dat[,groupingName]),n,replace=FALSE)
  stepsize = rep(floor(length(subj)/v),v)
  remainder = n-sum(stepsize)
  add = sample(1:v,remainder,replace=FALSE)
  stepsize = c(0,cumsum(stepsize + (1:v)%in%add))
  for(i in 1:v){
    testDat = dat[dat[,groupingName]%in%subj[(stepsize[i]+1):stepsize[i+1]],]
    data = dat[!dat[,groupingName]%in%subj[(stepsize[i]+1):stepsize[i+1]],]
    mod <- longRPart(model$nlmeFormula,model$rPartFormula,model$randomFormula,data,R=model$R,control=model$control)
    pred = predict(mod,testDat)
    # need to get the predicted lines at each node,
    # NOT just the regression coefficient
    terms = attr(terms(mod$nlmeFormula),"term.labels")
    timeVar = data[,names(data)==terms[1]]
    responseName = attr(terms(getResponseFormula(mod$lmeFormula)),"term.labels")
    continuous = !is.factor(timeVar)
    nodes = unique(mod$where)
    timeValues = unique(timeVar)
    for(j in 1:length(nodes)){
      subData = data[mod$where==nodes[j],]
      subTest = testDat[pred==mod$frame$yval[nodes[j]],]
      mod2 <- lme(mod$nlmeFormula,data=subData,random=mod$randomFormula,correlation=mod$R,na.action=na.omit)
      responseName = attr(terms(getResponseFormula(model$nlmeFormula)),"term.labels")
      for(k in unique(subTest[,names(subTest)==mod$groupingName])){
        observed <- subTest[subTest[,names(subTest)==mod$groupingName]==k,]
        # creating the fake patient (PROBLEMS HERE)
        # Things needed to be fixed:
        # 1. deal with non-standard time intervals
        # 2. deal with confounding
        expected <- observed
        expected[,names(expected)==mod$groupingName] <- -1
        subTimeVar <- unique(sort(mod2$data[!is.na(mod2$data[,names(mod2$data)==responseName]),names(mod2$data)==terms[1]]))
        expected[expected[,names(expected)==terms[[1]]]%in%subTimeVar,names(expected)==responseName] = as.vector(mod2$coeff$fixed)
        mod3 = lme(mod$nlmeFormula,data=rbind(observed,expected),random=mod$randomFormula,correlation=mod$R,na.action=na.omit)
        deviances[z,] <- c(k,-2*mod3$logLik)
        z = z+1
      }
    }
  }
  deviances
}

# Calculate a confidence interval for the first split using bootstrapping
#
#
lrpCI <- function(model,B=10,alpha=0.05){
  # idea is to take B bootstrap samples of size n, getting an optimal spliting
  # value at each node.  The CI makes sense for the first split, but for the
  # rest of the splits it is not as intuitive.  !!!For now, we'll just use the
  # first split, and discuss other options later
  R = model$R
  evaluation = model$functions$eval
  split = model$functions$split
  initialize = model$functions$init
  models = list()
  splitVal = rep(0,B)
  newControl = model$control
  newControl$maxdepth = 1
  newControl$cp=0
  n=length(unique(model$data[,model$groupingName]))
  for(b in 1:B){
    ind = sample(1:n,n,replace=TRUE)
    newDat = model$data[model$data[,model$groupingName]==ind[1],]
    newDat[,names(newDat)==model$groupingName] = rep(1,dim(newDat)[1])
    for(i in 2:n){
      tempDat = model$data[model$data[,model$groupingName]==ind[i],]
      tempDat[,names(newDat)==model$groupingName] = rep(i,dim(tempDat)[1])
      newDat = rbind(newDat,tempDat)
    }
    newDat[,names(newDat)==model$groupingName]=as.factor(newDat[,names(newDat)==model$groupingName])
    newModel = rpart(paste(model$groupingName,"~",row.names(model$splits)[1]),
                     data=newDat,method=list(eval=evaluation,split=split,
                     init=initialize),parms=model$data,control=newControl)
    models[[b]]=newModel
    if(!is.null(newModel$splits)){
      splitVal[b] = newModel$splits[1,4]
    }
    else{
      splitVal[b] = NA
    }
  }
  return(list(CI=quantile(splitVal,probs=seq(alpha/2,1-alpha/2,length.out=2),na.rm=TRUE),splitVal))
}

#
#
# Calculate the p-value using a permutation test
#
#
lrpPVal <- function(model,J){
  # idea is to permute the predictor while holding the response
  # fixed.  Do this J times, and then the proportion of improvements in deviance
  # that are better than the observed improvement is the p-val
  R = model$R
  Dobs = model$frame$dev[1] - (model$frame$dev[row.names(model$frame)==2]+model$frame$dev[row.names(model$frame)==3])
  evaluation = model$functions$eval
  split = model$functions$split
  initialize = model$functions$init
  models = list()
  D = rep(0,J)
  DFirst = rep(0,J)
  sub = length(unique(model$data[,model$groupingName]))
  newControl = model$control
  newControl$maxdepth=1
  for(j in 1:J){
    newDat = model$data
    ind = match(match(newDat[,model$groupingName],sample(1:sub,sub,replace=FALSE)),newDat[,model$groupingName])
    newDat[,names(newDat)==row.names(model$splits)[1]] = newDat[ind,names(newDat)==row.names(model$splits)[1]]
    newModel = rpart(paste(model$groupingName,"~",row.names(model$splits)[1]),data=newDat,method=list(eval=evaluation,split=split,init=initialize),parms=model$data,control=newControl)
    D[j] = newModel$frame$dev[1] - newModel$frame$dev[2] - newModel$frame$dev[3]
    models[[j]] = newModel
  }
  return(list(Dobs,D))
}

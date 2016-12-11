#library(rpart)
#library(nlme)
#
#
# Tree Building Function
#
#
longRPart <- function(lmeFormula,rPartFormula,randomFormula,data,weight=NULL,R=NULL,control = rpart.control()){
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
    split <- function(y, wt, x, parms, continuous){
      print(paste("splitting:",length(unique(x)),"values"))
      dev = vector()
      xUnique = unique(x)
      rootDev = lme(lmeFormula,data=parms[groupingFactor%in%y,],random=randomFormula,correlation=R,na.action=na.omit)$logLik
      #
      # for continuous variables
      # 
      if(continuous){
        # need to find splits for the UNIQUE values of x
        # no point in recalculating splits for the
        # same age/height/etc...
        for(i in xUnique){
          yLeft = y[x<=i]
          yRight = y[x>i]
          # build the models only if the split created satisfies the minbucket.
          # could add other split controls later
          if(length(yLeft) < control$minbucket || length(yRight) < control$minbucket){
            dev = c(dev,0)
          }
          else{
            modelLeft = try(lme(lmeFormula,data=parms[groupingFactor%in%yLeft,],random=randomFormula,correlation=R,na.action=na.omit),silent=TRUE)
            modelRight = try(lme(lmeFormula,data=parms[groupingFactor%in%yRight,],random=randomFormula,correlation=R,na.action=na.omit),silent=TRUE)
            if(class(modelLeft)=='lme' && class(modelRight)=='lme'){
              dev = c(dev,modelLeft$logLik+modelRight$logLik)
            }
            else{
              dev = c(dev,0)
            }
          }
        }
        # need to duplicate the observations for duplicate values of x
        good = rep(0,length(x))
        for(i in 1:length(xUnique)){
          good[x==xUnique[i]]=dev[i]
        }
        good = good[1:(length(good)-1)]
        list(goodness=good+abs(rootDev)*(good!=0)*2,direction=rep(-1,length(good)))
      }
      #
      ###for categorical variables
      #
      else{
        order = rep(0,length(xUnique))
        response = parms[,names(parms)==responseName]
        # establishing the ordering
        for(i in 1:length(xUnique)){
          order[i] = mean(response[x==xUnique[i]],na.rm=TRUE)
        }
        dir = sort(order,index.return=TRUE)$ix
        # testing the direction
        for(i in 1:(length(dir)-1)){
          yLeft = y[x%in%dir[1:i]]
          yRight = y[x%in%dir[(i+1):length(dir)]]
          # build the models only if the split created satisfies the minbucket.
          # could add other split controls later
          if(length(yLeft) < control$minbucket || length(yRight) < control$minbucket){
            dev = c(dev,0)
          }
          else{
            modelLeft = try(lme(lmeFormula,data=parms[groupingFactor%in%yLeft,],random=randomFormula,correlation=R,na.action=na.omit),silent=TRUE)
            modelRight = try(lme(lmeFormula,data=parms[groupingFactor%in%yRight,],random=randomFormula,correlation=R,na.action=na.omit),silent=TRUE)
            if(class(modelLeft)=='lme' && class(modelRight)=='lme'){
              dev = c(dev,modelLeft$logLik+modelRight$logLik)
            }
            else{
              dev = c(dev,0)
            }
          }
        }
        list(goodness=dev+abs(rootDev)*(dev!=0)*2,direction=dir)
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
    
    model = rpart(paste(groupingName,c(rPartFormula)),method=list(eval=evaluation,split=split,init=initialize),control=control,data=data,parms=data)
    model$lmeModel = lme(lmeFormula,data=data,random=randomFormula,correlation=R,na.action=na.omit)
    model$lmeFormula = lmeFormula
    model$randomFormula = randomFormula
    model$R = R
    model$data = data
    model$groupingName = groupingName
    model$rPartFormula = rPartFormula
    return(model)
  }
    

#' Plot the tree
#'
#' @param model Model object from longRPart2()
#' @param use.n Print N in each node
#' @param colors Colors to use
#' @param place Where to place the plot
#' @export

lrpTreePlot <- function(model,use.n=TRUE,colors=NULL,place="bottomright"){
  indexes = plot(model)
  confounding = (length(model$lmeModel$contrasts)>1)
  n = length(model$nodeLines)
  t = length(model$levels)
  plot(model,ylim=range(indexes$y)-c(diff(range(indexes$y))*.125,0),xlim=range(indexes$x)+c(-1,1)*diff(range(indexes$x)*.1))
  f = model$frame
  leaves = (1:dim(f)[1])[f$var=="<leaf>"]
  terms = attr(terms(model$lmeFormula),"term.labels")
  timeVar = model$data[,names(model$data)==terms[1]]
  responseName = attr(terms(getResponseFormula(model$lmeFormula)),"term.labels")
  continuous = !is.factor(timeVar)
  nodes = unique(model$where)
  timeValues = unique(timeVar)
  # need to come up with the plotting range here
  # plot(0,type='n',xlab='',ylab='',main='',xlim=c(1,8),ylim=c(0,110))
  plotList = list()
  for(i in 1:length(nodes)){
    dat=model$data[model$where==nodes[i],]
    mod = lme(model$lmeFormula,data=dat,random=model$randomFormula,correlation=model$R,na.action=na.omit)
    if(length(terms)>1){
      form = formula(paste(responseName,'~as.numeric(',terms[1],')|',paste(terms[-1],sep='\ '),sep=''))
    }
    else{
      form = formula(paste(responseName,'~as.numeric(',terms[1],')',sep=''))
    }
    plotDat = plot(mod,form)
    condLevels = plotDat$condlevels
    pts = list()
    for(j in 1:length(plotDat$panel.args)){
      pts[[j]] = tapply(plotDat$panel.args[[j]]$y,plotDat$panel.args[[j]]$x,mean)
    }
    plotList[[i]]=pts
  }
  subYLim = range(plotList,na.rm=TRUE)
  xRange = range(indexes$x)
  yRange = range(indexes$y)
  if(is.null(colors)){
    colors = rainbow(dim(f)[1],v=0.9)
  }
  for(i in 1:length(leaves)){
    limx = indexes$x[leaves[i]]+c(-1,1)*diff(xRange)*.075
    limy = indexes$y[leaves[i]]+c(-2,0)*diff(yRange)*.075
    rect(limx[1],limy[1],limx[2],limy[2],density=-1,col=gray(0.9))
    for(j in 1:length(plotList[[i]])){
      pts = plotList[[i]][[j]]
      pts = (pts - subYLim[1])/(subYLim[2]-subYLim[1])
      pts = pts*(limy[2]-limy[1])+limy[1]
      points(x=seq(limx[1],limx[2],length.out=length(pts)),pts,lwd=2,type='l',lty=j,col=colors[i])
    }
  }
  text(model,use.n=use.n)
  lineTypes = c("dashed", "dotted", "dotdash", "longdash", "twodash")
  text2 = character()
  if(length(terms)>1){
    for(i in 1:length(condLevels)){
      text2=c(text2,paste(lineTypes[i],"line =",condLevels[[i]][2]))
    }
    legend(place,legend=text2,pch=16,col='black',bty='n')
  }


}

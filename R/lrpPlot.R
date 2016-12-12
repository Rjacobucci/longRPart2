#' @export

lrpPlot <- function(model,smoothing="n",color=NULL,place="bottomright"){
  #smoothing can be either n for none, p for partial, or c for complete
  ### ONLY 'n' implemented for now
  terms = attr(terms(model$lmeFormula),"term.labels")
  timeVar = model$data[,names(model$data)==terms[1]]
  responseName = attr(terms(getResponseFormula(model$lmeFormula)),"term.labels")
  continuous = !is.factor(timeVar)
  nodes = unique(model$where)
  timeValues = unique(timeVar)
  if(is.null(color)){
    color = rainbow(length(nodes),v=0.8)
  }
  # need to come up with the plotting range here
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
    print(condLevels)
    pts = list()
    for(j in 1:length(plotDat$panel.args)){
      pts[[j]] = tapply(plotDat$panel.args[[j]]$y,plotDat$panel.args[[j]]$x,mean)
    }
    plotList[[i]]=pts
  }
  plot(as.numeric(timeValues),1:length(timeValues),type='n',xlab='',ylab='',main='',ylim=range(plotList))
  for(i in 1:length(plotList)){
    for(j in 1:length(plotList[[i]])){
      points(as.numeric(names(plotList[[i]][[j]])),plotList[[i]][[j]],type='l',lty=j,col=color[i])
    }
  }
  text1 = paste("node",nodes)
  lineTypes = c("dashed", "dotted", "dotdash", "longdash", "twodash")
  text2 = character()
  if(length(terms)>1){
    for(i in 1:length(condLevels)){
      text2=c(text2,paste(lineTypes[i],"line =",condLevels[[i]][2]))
    }
  }
  legend(place,legend=c(text1,text2),pch=16,col=c(color,rep('black',length(terms[-1]))))
}

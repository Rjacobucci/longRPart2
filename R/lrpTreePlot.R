#' Plot the tree
#'
#' @param model Model object from longRPart2()
#' @param use.n Print N in each node
#' @param colors Colors to use
#' @param place Where to place the plot
#' @export

lrpTreePlot <- function(model,use.n=TRUE,colors=NULL,place="bottomright"){


  if(model$method == "nlme"){
    indexes = plot(model$rpart_out)
    confounding = (length(model$lmeModel$contrasts) > 1)
    n = length(model$rpart_out$nodeLines)
    t = length(model$rpart_out$levels)
    plot(model$rpart_out, ylim = range(indexes$y) - c(diff(range(indexes$y)) *
                                                        0.125, 0), xlim = range(indexes$x) + c(-1, 1) * diff(range(indexes$x) *
                                                                                                               0.1))
    f = model$rpart_out$frame
    leaves = (1:dim(f)[1])[f$var == "<leaf>"]

    terms = attr(terms(model$nlme.model), "term.labels")
    responseName = attr(terms(getResponseFormula(model$nlme.model)),
                        "term.labels")
    timeVar = model$data[, names(model$data) == terms[1]]
    continuous = !is.factor(timeVar)
    nodes = unique(model$rpart_out$where)
    timeValues = unique(timeVar)
    plotList = list()



    for (i in 1:length(nodes)) {
      dat = model$data[model$rpart_out$where == nodes[i], ]


      mod = nls(formula = model$nlme.model, data = dat, start = model$start)

      ####

      for(k in 2:length(all.vars(model$nlme.model))){
        for(m in ncol(dat):1){
          if(all.vars(model$nlme.model)[k] == colnames(dat)[m]){
            time.metric = (all.vars(model$nlme.model)[k])
          }
        }
      }

      r1 = min(model$data[[time.metric]], na.rm = TRUE)
      r2 = max(model$data[[time.metric]], na.rm = TRUE)

      x = seq(from = r1, to = r2, length.out = 1000)


      new.df <- data.frame(time.metric,x)
      colnames(new.df) = c("NON",paste(new.df[1,1]))
      new.df[,"NON"] = NULL

      pts = list()

      pts[[1]] = tapply(predict(mod, new.df), new.df,
                        mean)

      plotList[[i]] = pts
    }
    plot(model$rpart_out, ylim = range(indexes$y) - c(diff(range(indexes$y)) *
                                                        0.125, 0), xlim = range(indexes$x) + c(-1, 1) * diff(range(indexes$x) *
                                                                                                               0.1))

    subYLim = range(plotList, na.rm = TRUE)
    xRange = range(indexes$x)
    yRange = range(indexes$y)
    if (is.null(colors)) {
      colors = rainbow(dim(f)[1], v = 0.9)
    }
    for (i in 1:length(leaves)) {
      limx = indexes$x[leaves[i]] + c(-1, 1) * diff(xRange) *
        0.075
      limy = indexes$y[leaves[i]] + c(-2, 0) * diff(yRange) *
        0.075
      rect(limx[1], limy[1], limx[2], limy[2], density = -1,
           col = gray(0.9))
      for (j in 1:length(plotList[[i]])) {
        pts = plotList[[i]][[j]]
        pts = (pts - subYLim[1])/(subYLim[2] - subYLim[1])
        pts = pts * (limy[2] - limy[1]) + limy[1]
        points(x = seq(limx[1], limx[2], length.out = length(pts)),
               pts, lwd = 2, type = "l", lty = j, col = colors[i])
      }
    }
    text(model$rpart_out, use.n = use.n)
    lineTypes = c("dashed", "dotted", "dotdash", "longdash",
                  "twodash")
    text2 = character()


  }

  else{

    indexes = plot(model$rpart_out)
    confounding = (length(model$lmeModel$contrasts) > 1)
    n = length(model$rpart_out$nodeLines)
    t = length(model$rpart_out$levels)
    plot(model$rpart_out, ylim = range(indexes$y) - c(diff(range(indexes$y)) *
                                                        0.125, 0), xlim = range(indexes$x) + c(-1, 1) * diff(range(indexes$x) *
                                                                                                               0.1))
    f = model$rpart_out$frame
    leaves = (1:dim(f)[1])[f$var == "<leaf>"]

    terms = attr(terms(model$fixedFormula), "term.labels")
    responseName = attr(terms(getResponseFormula(model$fixedFormula)),
                        "term.labels")

    timeVar = model$data[, names(model$data) == terms[1]]
    continuous = !is.factor(timeVar)
    nodes = unique(model$rpart_out$where)
    timeValues = unique(timeVar)
    plotList = list()




    for (i in 1:length(nodes)) {
      dat = model$data[model$rpart_out$where == nodes[i], ]


      mod = lm(model$fixedFormula, data = dat)

      if (length(terms) > 1) {
        form = formula(paste(responseName, "~as.numeric(",
                             terms[1], ")|", paste(terms[-1], sep = " "),
                             sep = ""))
      }
      else {
        form = formula(paste(responseName, "~as.numeric(",
                             terms[1], ")", sep = ""))
      }



      r1 = min(model$data[[terms[1]]], na.rm = TRUE)
      r2 = max(model$data[[terms[1]]], na.rm = TRUE)

      x = seq(from = r1, to = r2, length.out = 1000)


      new.df <- data.frame(terms[1],x)
      colnames(new.df) = c("NON",paste(new.df[1,1]))
      new.df[,"NON"] = NULL

      pts = list()

      pts[[1]] = tapply(predict(mod, new.df), new.df,
                        mean)

      plotList[[i]] = pts

    }

    subYLim = range(plotList, na.rm = TRUE)
    xRange = range(indexes$x)
    yRange = range(indexes$y)
    if (is.null(colors)) {
      colors = rainbow(dim(f)[1], v = 0.9)
    }
    for (i in 1:length(leaves)) {
      limx = indexes$x[leaves[i]] + c(-1, 1) * diff(xRange) *
        0.075
      limy = indexes$y[leaves[i]] + c(-2, 0) * diff(yRange) *
        0.075
      rect(limx[1], limy[1], limx[2], limy[2], density = -1,
           col = gray(0.9))
      for (j in 1:length(plotList[[i]])) {
        pts = plotList[[i]][[j]]
        pts = (pts - subYLim[1])/(subYLim[2] - subYLim[1])
        pts = pts * (limy[2] - limy[1]) + limy[1]
        points(x = seq(limx[1], limx[2], length.out = length(pts)),
               pts, lwd = 2, type = "l", lty = j, col = colors[i])
      }
    }
    text(model$rpart_out, use.n = use.n)
    lineTypes = c("dashed", "dotted", "dotdash", "longdash",
                  "twodash")
    text2 = character()
    if (length(terms) > 1) {
      for (i in 1:length(condLevels)) {
        text2 = c(text2, paste(lineTypes[i], "line =", condLevels[[i]][2]))
      }
      legend(place, legend = text2, pch = 16, col = "black",
             bty = "n")
    }



  }

}

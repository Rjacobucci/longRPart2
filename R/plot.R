#' Plot function for longRPart2
#'
#' @param x A model from lrp.
#' @param box.palette Color scheme for rpart.plot
#' @param ... Other arguments.
#' @method plot lrp
#' @export

plot.lrp <- function (x,box.palette="auto",...){

  rp = x$rpart_out

  node.fun1 <- function(x, labs, digits, varlen)
  {
    paste("dev", round(summary(x)$frame$dev,2))
  }

  rpart.plot::rpart.plot(rp,nn=T,box.palette=box.palette,node.fun=node.fun1)
}

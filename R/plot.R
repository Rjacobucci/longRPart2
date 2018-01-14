#' Plot function for longRPart2
#'
#' @param x A model from lrp.
#' @param ... Other arguments.
#' @method plot lrp
#' @export

plot.lrp <- function (x, ...){

  rp = x$rpart_out

  node.fun1 <- function(x, labs, digits, varlen)
  {
    paste("dev", round(summary(x)$frame$dev,2))
  }

  rpart.plot::rpart.plot(rp,nn=T,node.fun=node.fun1)
}

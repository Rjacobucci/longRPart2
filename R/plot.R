#' Plot function for longRPart2
#'
#' @param x A model from lrp.
#' @param ... Other arguments.
#' @method plot lrp
#' @export

plot.lrp <- function (x, ...){
  plot(x$rpart_out);text(x$rpart_out)
}

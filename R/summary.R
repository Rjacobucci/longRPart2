#' Summary results from lrp.
#'
#' @param object An object from lrp.
#' @param ... Other arguments.
#' @method summary lrp
#' @export
#'


summary.lrp <- function(object, ...){

  ret <- list()

  ret$importance <- object$rpart_out$variable.importance


  mat <- list()
  count=1
  for(i in 1:length(object$fixed_effects)){
    if(is.null(object$fixed_effects[[i]])){
      count=count
    }else{
      if(count==1){
        names <- c(names(object$fixed_effects[[i]]))
      }
      mat[[count]] <- c(i,object$fixed_effects[[i]],object$resid.var[[i]])
      count = count+1
    }
  }

  mm = data.frame(matrix(unlist(mat),nrow=(count-1),byrow=T))
  colnames(mm)[1] <- "node"
  names(mm)[2:(ncol(mm)-1)] <- names
  names(mm)[ncol(mm)] <- "resid"

  ret$parameters <- mm

  class(ret) <- "summary.lrp"
  print(ret)
}

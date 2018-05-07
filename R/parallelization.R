#' (Internal function) Check the different value of X and Y, internally used for
#' parallelization
#'
#' @param x The first number
#' @param y The second number
#'
#' @return The different number of \code{x} and \code{y}.
#'
#' @examples
#' print(diff.xy(2.4,6.7))
#' print(diff.xy(8.23,3.5))
#' print(diff.xy(4.3,4.3))
diff.xy = function(x,y){
  ret = abs(x-y)
  return(ret)
}

#' (Internal function) Calculae a vector of EigenFit values, internally used for
#' parallelization
#'
#' @param eigen.value A vector of Eigenvalues return from \code{\link{svd}}
#' ($d), \code{\link{rARPACK::svds}} ($d),  \code{\link{eigen}}
#' ($values) or \code{\link{rARPACK::eigs}} ($values).
#'
#' @return A vector of all possible values for EigenFit.
#'
#' @examples
#'  X <- sort(runif(10, min = 0, max = 3), decreasing = TRUE)
#'  print(cal.eigen.fit(X))
cal.eigen.fit = function(eigen.value){
  I=log(eigen.value)
  X = I
  Y = I
  Y = Y[2:length(I)]
  X = X[1:(length(I)-1)]
  ret = mapply(diff.xy,X,Y)
  return(ret)
}

#' (Internal function) Calculate a vector of different values from a vector of
#' EigenFit values, internally used for parallelization
#'
#' @param eigen.value A vector of Eigenvalues return from \code{\link{svd}}
#' ($d), \code{\link{rARPACK::svds}} ($d),  \code{\link{eigen}}
#' ($values) or \code{\link{rARPACK::eigs}} ($values).
#'
#' @return A vector of different values from a vector of EigenFit values
#'
#' @examples
#'  X <- sort(runif(10, min = 0, max = 3), decreasing = TRUE)
#'  print(diff.eigen.fit(X))
diff.eigen.fit = function(eigen.value){
  I=log(eigen.value)
  deviation=I-median(I)
  I=I[which(deviation>0)]

  X = I
  Y = I
  Y = Y[2:length(I)]
  X = X[1:(length(I)-1)]
  ret = log(mapply(diff.xy,X,Y))
  return(ret)
}

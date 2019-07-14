#library(tidyverse)
#library(dplyr)

#' @title Create Cluster blal
#' @param y a  (nx2) matrix (nx2) First column: number of bin, Second column: number of observation per bin
#' @param k a numeric value - number of components (fitted distributions)
#' @param method the method to be used for finding the k groups of data. For more information see Details

#' @return kx2 matrix with mean and variance per component
#' @description The function estimates the mean and variance of the k groups of the underlying inputdata y.
#' @details For the input parameter method are two options available, namely binbased and quantiles.
#' binbased: the group calculation is based on bins. As a result of that each group has the same number of bins
#' quantiles: the group calculation is based on the observed values per bin. As a result of that each group has the same number of observations.
#' @examples
#' v <- c(2, 4, 5,6,5,2,2, 1, 1, 2,  2, 1,6,7,8,7,6, 5, 2,1)
#'
#' data <- data.frame(name = 1:length(v)+5, v)
#' data  <- as.matrix(data)
#' createCluster(as.matrix(data), 2, method = 'binbased')

#' @export
createCluster <- function(y,k, method = "quantile"){

  #y matrix with columns (name of bin, number of observation)
  #k number of groups
  # method quantiles and binbased

 # if(class(y) != "matrix") warning("y is not a matrix")
 # if(class(k) != "numeric") warning("k is not a numeric vector")
  if(any(k<0)) warning("only positiv x values are allowed")
  if(any(y<0)) warning("only positiv y values are allowed")

  # Delete 0 Observations
  y <- y[-1, ]
  y <- y[y[,2] != 0, ]

  parameters <- data.frame(mu = c(1:k),sigma2 = c(1:k))


  if (method == 'binbased'){
    groups <- round(length(y[,1])/k)
    y <- data.frame(y)


    lowerBound <- 0
    upperBound <- groups

    # GRUEN: Assign group with rep
    # GRUEN: estimates for pi
    # GRUEN: Delete 0
    for(i in 1:k){

   #   yNew <- y %>% dplyr::filter(y[,1] <= upperBound & lowerBound < y[,1])
      yNew <- y[y[,1] <= upperBound & lowerBound < y[,1], ]

      parameters[i,1] <- mean(rep(yNew[,1],yNew[,2]))
      parameters[i,2] <- stats::var(rep(yNew[,1],yNew[,2]))

      lowerBound <- upperBound
      upperBound <- upperBound + groups

    }

  }

  if ( method == "quantile"){
    groups <- round(sum(y[,2])/k)

    y.splitted <- rep(y[,1],y[,2])
    g <- rep(1:k, each= groups)

    df <- data.frame(d = y.splitted, group = g[1:length(y.splitted)] )
    parameters[,1] <- stats::aggregate(df$d, list(df$group), FUN = 'mean')$x
    parameters[,2] <- stats::aggregate(df$d, list(df$group), FUN = 'var')$x

  }


  return(parameters)
}

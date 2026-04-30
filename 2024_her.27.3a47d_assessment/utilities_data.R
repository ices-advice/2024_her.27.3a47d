## Put data into FLR format

library(reshape2)  # dcast

formatLAI <- function(x,minYr,maxYr)
{
  x <- dcast(x[,c("Year","LAIUnit","L..9")], Year~LAIUnit, value.var="L..9")
  rownames(x) <- x$Year + 1900
  missingYears <- (minYr:maxYr)[which(!minYr:maxYr %in% rownames(x))]
  if(length(missingYears) > 0){
    xextra <- cbind(Year=missingYears,
                    do.call(cbind,
                            lapply(as.list(1:(ncol(x)-1)),
                                   function(y){
                                     return(rep(NA,length(missingYears)))
                                   })))
    rownames(xextra) <- missingYears
    colnames(xextra) <- colnames(x)
    x <- rbind(x,xextra)
  }
  x <- x[sort.int(rownames(x),index.return=TRUE)$ix,]
  x <- x[, -grep("Year",colnames(x))]
  colnames(x) <- 0:(ncol(x)-1)
  x <- as.matrix(x)
  attr(x,"time") <- c(0.67, 0.67)
  return(x)
}

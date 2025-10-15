#' Creates random samples and produces barplots
#'
#' @param n an integer that describes the size of each sample
#' @param iter an integer that describes the amount of times a trial should be run
#' @param time time between each iteration
#' @importFrom grDevices rainbow
#'
#' @returns
#' @export
#'
#' @examples
#' mysample(n=1000, iter=30, time = 1)
mysample=function(n, iter=10,time=0.5){
  for( i in 1:iter){
    #make a sample
    s=sample(1:10,n,replace=TRUE)
    # turn the sample into a factor
    sf=factor(s,levels=1:10)
    #make a barplot
    barplot(table(sf)/n,beside=TRUE,col=rainbow(10),
            main=paste("Example sample()", " iteration ", i, " n= ", n,sep="") ,
            ylim=c(0,0.2)
    )

    #release the table
    Sys.sleep(time)
  }
}

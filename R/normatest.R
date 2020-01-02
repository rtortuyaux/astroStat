#adapted from RQuery (sthda.com)

#' norma test
#'
#' @param x
#' @param graph
#' @param colname
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
normaTest <- function(x, graph = TRUE, colname = colname,...)
{
  # Significance test
  #++++++++++++++++++++++
  shapiro.p<-signif(shapiro.test(x)$p.value,1)

  if(graph){
    # Plot : Visual inspection
    #++++++++++++++++
    h<-hist(x, col="lightblue", main="Histogram",
            xlab="Data values", ...)
    m<-round(mean(x),1)
    s<-round(sd(x),1)
    mtext(paste0("Mean : ", m, "; SD : ", s),
          side=3, cex=0.8)
    # add normal curve
    xfit<-seq(min(x),max(x),length=40)
    yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
    yfit <- yfit*diff(h$mids[1:2])*length(x)
    lines(xfit, yfit, col="red", lwd=2)
    # qq plot
    qqnorm(x, pch=19, frame.plot=FALSE,main="Normal Q-Q Plot")
    qqline(x)
    mtext(paste0("Shapiro-Wilk, p-val : ", shapiro.p),
          side=3, cex=0.8)
    #create new folder with plots
    dev.print(device = pdf, file = paste0(getwd(), "/NormaTest/", "Norma_", colname, ".pdf"), width = 7, height=7,)
  } else {
    unlink("NormaTest/", recursive = TRUE)
  }
  return(shapiro.p)
}

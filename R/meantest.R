#' mean test
#'
#' @param x
#' @param y
#' @param paired
#' @param graph
#' @param colname
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

meanTest <- function(x, y = NULL, paired = FALSE, graph = TRUE, colname = colname, ...)
{

  #Preliminary test : normality and variance tests
  var.equal = FALSE
  cat("\n", "pre-test conditions :", "\n", sep = "")

  #One sample t test
  if(is.null(y)){
    if(graph) par(mfrow=c(1,2))
    shapiro.px<-normaTest(x, graph, colname=colname)
    if(shapiro.px < 0.05)
      cat(colnames(x), "is not normally distributed : Shapiro-Wilk test p-value = ", shapiro.px, "\n", sep = " ")
  }

  #Two samples t test
  if(!is.null(y)){

    #if unpaired t test
    if(!paired){
      if(graph) par(mfrow=c(2,2))
      #normality test
      shapiro.px<-normaTest(x, graph, colname=colname)
      shapiro.py<-normaTest(y, graph,colname=colname)
      if(shapiro.px < 0.05 | shapiro.py < 0.05){
        cat(colnames(x), "(x or y) is not normally distributed :",
            "Shapiro test p-value :", shapiro.px,
            "(for x) and", shapiro.py, "(for y)", "\n", sep = " ")
      }
      #Check for equality of variances
      if(var.test(x,y)$p.value >= 0.05) {
        var.equal=TRUE
      }
      if((var.test(x,y)$p.value < 0.05) & ((shapiro.px >= 0.05 & shapiro.py >= 0.05) | (NROW(x) > 30) & (NROW(y) > 30))) {
        cat("Variance of", colnames(x), "are different : use Welch-Satterthwaite test", "\n", sep = " ")
      }
      if((var.test(x,y)$p.value >= 0.05) & ((shapiro.px >= 0.05 & shapiro.py >= 0.05) | (NROW(x) > 30) & (NROW(y) > 30))) {
        cat("Variance of", colnames(x), "are not statistically different : use Student t-test", "\n", sep = " ")
      }
    }

    #Paired t-test
    else {
      if(graph) par(mfrow=c(1,2))
      d = x-y
      shapiro.pd<-normaTest(d, graph,colname=colname)
      if(shapiro.pd < 0.05 )
        cat("The difference d (x-y) for", colnames(x), "is not normally distributed :",
            " Shapiro-Wilk test p-value : ", shapiro.pd, "\n", sep = " ")
    }
  }

  #central limit theorem
  if ((NROW(x) > 30) & (NROW(y) > 30)) {
    cat("n > 30 in each group : application of central limit theorem", "\n")
  }

  #Student's t-test
  if (paired) {
    if((shapiro.pd >= 0.05) | (NROW(x) > 30)) {
      res <- t.test(x, y, paired=paired, var.equal=var.equal, ...)
      return(res)}
    else {
      res <- wilcox.test(x, y, paired=paired, ...)
      return(res)
    }
  } else {
    if((shapiro.px >= 0.05 & shapiro.py >= 0.05) | (NROW(x) > 30) & (NROW(y) > 30)) {
      res <- t.test(x, y, paired=paired, var.equal=var.equal, ...)
      return(res) }
    else {
      res <- wilcox.test(x, y, paired=paired, ...)
      return(res)
    }
  }
}

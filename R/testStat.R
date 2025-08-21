#' calculate a test statistc
#'
#' For a given randomisation sequence and response, calculate a certain test statistic.
#' For now, we only look at the normalized difference of ranks
#'
#' @param R randomisation sequence.
#' @param y observed response.
#' @param type character string indicating the type of test statistic for the comparison. Can be
#' \describe{
#' \item{\code{"cdr"}}{for centralized difference in ranks (linear rank test),}
#' \item{\code{"lrs"}}{for linear rank test when the responses are already ranked,}
#' \item{\code{"lrs2"}}{for linear rank test when the responses are not yet ranked,}
#' \item{\code{"ds"}}{for difference in sums,}
#' \item{\code{dm}}{for difference in means,}
#' \item{\code{"se"}}{for sum of group E,}
#' \item{\code{"sc"}}{for sum of group C,}
#' \item{\code{"dmed"}}{for difference in medians, or}
#' \item{\code{"t"}}{for the t-test statistic.}
#'}
#'
#'
#'
#' @export
testStat<-function(R,y,type="lrs"){

  if(type=="lrs2"){ # linear rank statistic
    sum(R*cf(y))
  }

  else if(type=="lrs"){ # if the responses are already ranked
    # Temporary added for simulation purposes
    R <- R[!(y==0)]
    y <- y[!(y==0)]
    m <- mean(y, na.rm = T)
    sum(R*(y-m), na.rm = T)
  }

  else if(type=="cdr") { # difference of ranks
    sum((2*R-1)*cf(y))
  }

  else if (type=="ds"){
    2*sum(R*y) - sum(y)
  }

  else if(type=="dm") { # difference in means
    if(sum(R)%in%c(0,length(R)))
      0
    else
      sum(y[R==1])/sum(R==1)-sum(y[R==0])/sum(R==0)
  }

  else if(type=="se") { # sum of experimental group
    if(sum(R)==0)
      0
    else
      sum(y[R==1], na.rm=TRUE)
  }

  else if(type=="sc") { # sum of control group
    sum((1-R)*y)
  }

  else if(type=="dmed") {
    median(y[R==1])-median(y[R==0])
  }

  else if(type=="t") {
    n <- sum(R==1)
    m <- sum(R==0)
    x1 <- sum(y[R==1])/n
    x2 <- sum(y[R==0])/m
    s1sq <- var(y[R==1])
    s2sq <- var(y[R==0])

    sqrt(n*m/(n+m)) * (x1-x2)/sqrt(((n-1)*s1sq + (m-1)*s2sq)/(n+m-2))
  }

  else warning("Type of test statistic not supported")
}

#' Coefficients for the linear rank test
#'
#' @param y observed response
#'
cf <- function(y) {
  a <- rank(y)
  a-length(a)/2-0.5
}

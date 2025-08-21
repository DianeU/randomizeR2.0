#' @include util_condRefSet_mv.R

# --------------------------------------------
# Function for validity check
# --------------------------------------------

validateRandTest <- function(object) {
  errors <- character()
  lengthType <- length(object@type)
  if (lengthType != 1) {
    msg <-
      paste("Type is length ", lengthType, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  type <- object@type[1]
  #if (!(type %in% c("lrs","cdr", "dm", "se", "sc", "dmed", "t"))) {
  if (!(type %in% c("lrs2", "lrs","cdr", "ds", "se", "dm"))) {
    msg <-
      paste(
        #"(First) Argument of type is ", type, ". Should be \"lrs\", \"cdr\", \"dm\", \"se\", \"sc\", \"dmed\" or \"t\"."
        "(First) Argument of type is ", type, ". Should be \"lrs2\", \"lrs\", \"cdr\", \"ds\", \"dm\ ,\"se\"."
        , sep = ""
      )
    errors <- c(errors, msg)
  }
  if (length(errors) == 0)
    TRUE
  else
    errors
}


# --------------------------------------------
# Class definition for randomizationTest
# --------------------------------------------
#' Randomization Tests
#' 
#' Randomization tests yield the only valid test for the presence of treatment 
#' effects.
#' 
#' @slot pValue Probability of rejection
#' @slot obs observation consisting of response and randomization sequence
#' @slot type type of test statistic that is used
#' 
#' @references Rosenberger Lachin 2016
setClass("randomizationTest", slots = c(pValue = "numeric", obs = "numeric", 
                                        type = "character", K = "numeric"),
         validity = validateRandTest
         )


# --------------------------------------------
# Show function for randomizationTest
# --------------------------------------------
setMethod("show", "randomizationTest", function(object) {
  validObject(object)
  cat("\nObject of class \"", class(object)[1], "\"\n\n", sep = "")
  cat("         p-Value =", object@pValue, "\n")
  cat("         Observed test statistic =", object@obs, "\n")
  cat("         Type of test statistic =", object@type, "\n")
  cat("         Length of the reference set = ", object@K, "\n")
  cat("\n")
})



# ------------------------------------------------------
# Methods and Generic Functions for Randomization Tests
# ------------------------------------------------------

#' @title Conduct a Randomization Test
#'
#' @description Calculates the observed value of the test statistic and compares
#' it to the values it would have attained if one of the sequences of the
#' reference set had been observed. 
#'
#' @param obs Observations of the trial containing a randomization sequence and the responses, see \code{\link{genObs}}
#' @param refSet Reference set for the present comparison
#' @param type character string indicating the type of test statistic for the comparison. Can be 
#' \describe{
#' \item{\code{"lrs2"}}{linear rank test statistic,}
#' \item{\code{"cdr"}}{for centrlized difference in ranks (linear rank test),}
#' \item{\code{dm}}{for difference in means,}
#' \item{\code{"se"}}{for sum of group E,}
#' \item{\code{"sc"}}{for sum of group C,} 
#' \item{\code{"dmed"}}{for difference in medians, or}
#' \item{\code{"t"}}{for the t-test statistic.}
#'}
#' @examples
#' # from randomizeR
#' params <- rarPar(10)
#' rs <- genSeq(params)
#' endp <- normEndp(c(0,1), c(1,1))
#' refSet <- genSeq(params, 100)
#' # new in infeR
#' obs <- genObs(rs, endp)
#' randTest(obs, refSet)
#'
#' @return Returns the type of test statistic that had been used, the observed value of the test statistic and the the p-value.
#' @name conductRandTest
NULL


#' @rdname conductRandTest
#'
setGeneric("randTest", function(obs, refSet, type) standardGeneric("randTest"))


#' @rdname conductRandTest
#'
#' @export
setMethod("randTest", signature(obs = "observation", refSet = "rRandSeq", type = "character"),
          function(obs, refSet, type) {
            if (!(type %in% c("lrs2","lrs","cdr","ds","se", "dm")))
              return (stop("Invalid type, please check documentaion"))
            y <- obs@resp
            if(any(is.na(y)))
              refSet@M <- condRefSet(obs, refSet)
            S_obs <- testStat(obs@rs@M, y, type = type)
            s <- apply(refSet@M, 1, function(x) {
              testStat(x, y, type = type)
            })
            s <- c(s, S_obs)
            if(length(s)==1){
              pValue <- 1
            } else {
              pValue = getPval_obs_t(S_obs, s, rep(1 / (nrow(refSet@M) + 1), nrow(refSet@M) + 1)) 
            }
            new("randomizationTest", pValue = pValue, obs = S_obs, type = type, K = nrow(refSet@M))
          })


#' @rdname conductRandTest
#' 
#' @export
setMethod("randTest", signature(obs = "observation", refSet = "rRandSeq", type = "missing"),
          function(obs, refSet, type) {
            type <- "lrs"
            y <- obs@resp
            if(any(is.na(y)))
              refSet@M <- condRefSet(obs, refSet)
            S_obs <- testStat(obs@rs@M, y, type = type)
            s <- apply(refSet@M, 1, function(x) {
              testStat(x, y, type)
            })
            s <- c(s, S_obs)
            pValue = getPval_obs_t(S_obs, s, rep(1 / (nrow(refSet@M) + 1), nrow(refSet@M) + 1))
            new("randomizationTest", pValue = pValue, obs = S_obs, type = type, K = nrow(refSet@M))
          })


#' @rdname conductRandTest
#' 
#' @export
setMethod("randTest", signature(obs = "observation", refSet = "randSeq", type = "character"),
          function(obs, refSet, type) {
            if (!(type %in% c("lrs2","lrs","cdr", "ds","se", "dm")))
              return (stop("Invalid type, please check documentaion"))
            y <- obs@resp
            if(any(is.na(y)))
              refSet@M <- condRefSet(obs, refSet)
            S_obs <- testStat(obs@rs@M, y, type = type)
            s <- apply(refSet@M, 1, function(x) {
              testStat(x, y, type)
            })
            p <- getProb(refSet)
            pValue <- getPval_obs_t(S_obs, s, p)
            new("randomizationTest", pValue = pValue, obs = S_obs, type = type, K = nrow(refSet@M))
          })


#' @rdname conductRandTest
#' 
#' @export
setMethod("randTest", signature(obs = "observation", refSet = "randSeq", type = "missing"),
          function(obs, refSet, type) {
            type <- "lrs"
            y <- obs@resp
            if(any(is.na(y)))
              refSet@M <- condRefSet(obs, refSet)
            S_obs <- testStat(obs@rs@M, y, type = type)
            s <- apply(refSet@M, 1, function(x) {
              testStat(x, y, type)
            })
            p <- getProb(refSet)
            pValue <- getPval_obs_t(S_obs, s, p)
            new("randomizationTest", pValue = pValue, obs = S_obs, type = type, K = nrow(refSet@M))
          })
# --------------------------------------------
# Function for validity check
# --------------------------------------------

validateSimRTClass <- function(object) {}

# --------------------------------------------
# Class definition for randomizationTest
# --------------------------------------------
#' Randomization Tests
#' 
#' Randomization tests yield the only valid test for the presence of treatment 
#' effects.
#' 
#' @slot res Proportion of tests that had a result less than the alpha level.
#' @slot randPar Representation of the randomization procedure that is being simulated.
#' @slot issue Issue for the comparison. Can be of class selBias, mVal, or \code{"missing"}.
#' @slot endp Endpoint of the responses that are generated
#' @slot L Size of the MC reference set that is generated for each tet
#' @slot r Number of tests that are conducted
#' @slot alpha alpha level of the test
#' @slot type type of test statistic that is used
#' 
#' @references Rosenberger Lachin 2016
#' 
#' @name simulationRT-class
#' 
#' @exportClass simulationRT
setClass("simulationRT", slots = c(res = "numeric", randPar = "randPar", 
                                    issue = "character", endp = "endpoint", L = "numeric",
                                    r = "numeric", alpha = "numeric", type = "character"),
         validity = validateSimRTClass
)


# --------------------------------------------
# Show function for randomizationTest
# --------------------------------------------
setMethod("show", "simulationRT", function(object) {
  validObject(object)
  cat("\nObject of class \"", class(object)[1], "\"\n\n", sep = "")
  cat("         Percent below alpha =", object@res, "\n\n")
  cat("         Randomization Procedure used: \n")
  print(object@randPar)
  cat("\n")
  cat("         Issue used =", object@issue, "\n\n")
  cat("         Length of the reference sets =", object@L, "\n\n")
  cat("         Number of p-values considered =", object@r, "\n\n")
  cat("         Alpha =", object@alpha, "\n\n")
  cat("         Type of test statistic =", object@type, "\n\n")
  cat("\n")
})

# --------------------------------------------
# Check Input and warn user
# --------------------------------------------
validateInput <- function(endp, L, r, alpha, type, cores){
  if(alpha > 1 || alpha < 0 || !is.numeric(L))
    return (stop("Invalid input 'alpha'. Should be a number between 0 and 1"))
  if(!is(endp, "endpoint"))
    return (stop("Invalid input 'endp'. Should be an endpoint object"))
  if(!(type %in% c("lrs2","lrs","cdr", "ds", "dm")))
    return (stop("(First) Argument of type is ", type, ". Should be \"lrs2\", \"lrs\", \"cdr\", \"ds\", \"dm\"."))
  if(cores < 1 || cores > 64)
    return (stop("You can't use ", cores, " cores. Should use at least 1 and not more than 64."))
  if(cores %% 1 != 0)
    return (stop("Cores must be an integer, not a float number"))
}

# --------------------------------------------
# Generic function for simRT
# --------------------------------------------

#' @title Simulate Randomization Test
#' 
#' @description Simulates repeated randomization tests and calculates
#' the proportion of tests that yielded a p-value below the significance level.
#' 
#' @param randPar a \code{\link[randomizeR]{randPar}} object indicating the randomization 
#' procedure that should be simulated
#' @param endp \code{\link[randomizeR]{normEndp}} object indicating the endpoint of the responses
#' @param issue The observations can be affected by selection bias 
#' (\code{\link[randomizeR]{issue}}) or missing values
#' @param L integer that denotes the length of the reference set
#' @param r integer that denotes how many repetitions should be conducted
#' @param alpha significance level of the tests
#' @param type Type of test statistic used
#' @param cores Number of threads created to make the computation. Default is 1 - no parallelisation
#' 
#' @examples 
#' randPar <- rarPar(8)
#' issue <- selBias("CS", 0.05, "sim", 0.05)
#' simRT(randPar = randPar, issue = issue)
#' 
#' #For different iteration steps with more/less sequences, specify r and L
#' L <- 100
#' r <- 200
#' simRT(randPar = randPar, issue = issue, L = L, r = r)
#' 
#' @return proportion of tests that yielded a p-value below the significance level
#' 
#' @seealso Randomization procedures \link[randomizeR]{randPar}.
#' @seealso Endpoint of the responses \link[randomizeR]{normEndp}.
#' @seealso Responses might be influenced by bias, see \link[randomizeR]{issue}.
#' @name simulateRandTest
NULL

#' @rdname simulateRandTest
#' @export
setGeneric("simRT", function(randPar, issue, endp = normEndp(c(0,1), c(1,1)), L = 500, r = 100, alpha = 0.05, type = "lrs", cores = 1) standardGeneric("simRT"))


#' @rdname simulateRandTest
#' 
setMethod("simRT", signature(randPar = "randPar", issue = "missing"), 
        function(randPar, issue, endp, L, r, alpha, type, cores){
          validateInput(endp, L, r, alpha, type, cores)
          if(cores == 1){
            p_list <- lapply(1:r, function(x){
              rs <- genSeq(randPar)
              obs <- genObs(rs, endp)
              refSet <- genSeq(randPar, L)
              randTest(obs, refSet, type)@pValue
            })
          } else {
            cl <- makeCluster(cores)
            clusterCall(cl, function() library(infeR))
            clusterExport(cl, list("randPar","endp", "r", "L", "type"), envir = environment())
            p_list <- parLapply(cl = cl, x = 1:r, fun = function(x){
              rs <- genSeq(randPar)
              obs <- genObs(rs, endp)
              refSet <- genSeq(randPar, L)
              randTest(obs, refSet, type)@pValue
            })
            stopCluster(cl)
          }
          res <- sum(p_list <= alpha)/r
          new("simulationRT", res = res, randPar = randPar, issue = "None", endp = endp, L = L, r = r, alpha = alpha, type = type)
        }
)

#' @rdname simulateRandTest
#' 
setMethod("simRT", signature(randPar = "randPar", issue = "issue"), 
        function(randPar, issue, endp, L, r, alpha, type, cores){
          validateInput(endp, L, r, alpha, type, cores)
          if(cores == 1){
            p_list <- lapply(1:r, function(x){
              rs <- genSeq(randPar)
              obs <- genObs(rs, endp, issue)
              refSet <- genSeq(randPar, L)
              randTest(obs, refSet, type)@pValue
            })
          }
          else {
            cl <- makeCluster(2)
            clusterCall(cl, function() library(infeR))
            clusterExport(cl, list("randPar","endp", "issue", "r", "L", "type"), envir = environment())
            p_list <- parLapply(cl = cl, x = 1:r, fun = function(x){
              rs <- genSeq(randPar)
              obs <- genObs(rs, endp, issue)
              refSet <- genSeq(randPar, L)
              randTest(obs, refSet, type)@pValue
            })
            stopCluster(cl)
          }
          res <- sum(p_list <= alpha)/r
          if(class(issue)[1] == "selBias")
            issueString <- "Selection Bias"
          else if(class(issue)[1] == "chronBias")
            issueString <- "Chronological Bias"
          new("simulationRT", res = res, randPar = randPar, issue = issueString, endp = endp, L = L, r = r, alpha = alpha, type = type)
        }
)

#' Read entries from a csv file and create an observation object. The file must
#' contain just two columns, where the first consists of the randomization sequence
#' and the second of the corresponding responses.
#' 
#' @export 
readObs <- function(){
  # read the csv file
  newObs <- read.csv(choose.files(caption = "Select the observation file",
                                  multi = FALSE, filters = cbind("Only csv (*.csv)", "*.csv")), 
                                  header = FALSE, sep = ",", quote = "\"",
                                  dec = ".", fill = TRUE, comment.char = "")
  # get randomization Procedure
  randProc <- createParam(N = nrow(newObs))
  
  # run checks for:
  # More columns than two in the csv file
  if(ncol(newObs)>2)
    return(stop("Too many columns. You just need 2"))
  # If strings in the csv file
  if(!is.numeric(newObs[,1]) || !is.numeric(newObs[,2]))
    return(stop("Only numerics allowed in the csv file"))
  # If columns have unequal length
  if(length(newObs[,1]) != length(newObs[,2]))
    return(stop("Input file has different Column sizes"))
  r_sequence <- newObs[,1]
  response <- newObs[,2]

  rs <- genSeq(randProc)
  rs@M <- t(r_sequence)  

  obs <- new("observation", rs = rs, resp = response)
  return(obs)
}



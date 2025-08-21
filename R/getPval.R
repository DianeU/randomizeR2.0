#' calculate a exact p-value of non-parametric test
#' 
#' The p-value is the probability of obaining a test statistic as extreme or 
#' more extreme than the observed test statistic. More extreme is interpreted as
#' having a more extreme absolute value than the absolute value of the observed
#' test statistic.
#' 
#' @param ind index of the sequence of M that was actually used for the trial
#' @param ts vector of the test statistics for the sequences in the rows of M
#' @param probs vector with the probabilities of the sequences contained in M.
#' 
#' @export 
getPval<- function(ind, ts, probs) {
  sum(probs[abs(ts)>=abs(ts[ind])])
}


#' calculate a exact p-value for given observed value of the test statistic
#' 
#' The p-value is the probability of obaining a test statistic as extreme or 
#' more extreme than the observed test statistic. More extreme is interpreted as
#' having a more extreme absolute value than the absolute value of the observed
#' test statistic.
#' 
#' @param t_obs observed value of the test statistic
#' @param ts vector of the test statistics for the sequences in the rows of M
#' @param probs vector with the probabilities of the sequences contained in M.
#' 
#' @export 
getPval_obs_t<- function(t_obs, ts, probs) {
  sum(probs[abs(ts)>=abs(t_obs)])
}

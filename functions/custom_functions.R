## Custom functions

expit <- function(log_odds) {
  return(1 / (1 + exp(-log_odds)))
}

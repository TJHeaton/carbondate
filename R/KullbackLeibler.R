.KLD <- function(P, Q) {
  # P = P/sum(P)
  # Q = Q/sum(Q)
  return (sum(P * log(P / Q)))
}

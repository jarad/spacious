# Computes initial values for matern covariance matrix
initial.theta = function(y, S, p.nugget=0.1, p.range=0.1) {
  sill = var(y)              # 
  nugget = p.nugget*sill
  partial.sill = sill-nugget

  # Calculate <p.range>th percentile of distance between 1000 random points
  n = length(y)
  points = sample(n, min(n,1000))
  d = dist(S[points,])
  range = quantile(d,p.range)
  
  smoothness = 0.5

  return(c(nugget, partial.sill, 1/range, smoothness))
}

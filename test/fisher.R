if (0) {
# get basic MLE estimate for Poisson distribution using Fisher scoring
set.seed(311)
n <- 100
lambda <- 5
y <- rpois(n, lambda)

tol <- 1e-4
iters <- 20
alphas <- 10

#"fs" <- function(sum_y, n, oalpha) { oalpha - n + exp(-oalpha) * sum_y }
"fs" <- function(ybar, oalpha) { oalpha - 1 + exp(-oalpha) * ybar }

for (i in 1:(iters-1)) {
	#alphas[i+1] <- fs(sum(y), n, alphas[i])
	alphas <- c(alphas,fs(mean(y), alphas[i]))
	if ( abs(alphas[i+1]-alphas[i]) < tol) {
		print(c(alphas[i+1],alphas[i]))
		break;
	}
}

print(alphas)
print(exp(alphas))

}

if (1) {
# get basic MLE estimate for linear regression using Fisher scoring

set.seed(311)
n <- 500
sigma <- 1    # assume sigma known
beta <- c(0, 5)
X <- cbind(1, rnorm(n))
Y <- rnorm(n, mean=X %*% beta, sd=sigma)

"fs" <- function(Y, X, beta) {
	sx <- sum(X)
	sx2 <- sum(X^2)
	sy <- sum(Y)
	syx <- sum(Y*X)
	N <- length(Y)

	FI <- matrix( c(N, sx, sx, sx2), nrow=2)
	score <- matrix(c(sy - N * beta[1] - beta[2]*sx, syx - beta[1]*sx - beta[2] * sx2), nrow=2)
#print(FI)
#print(score)

	beta + chol2inv(chol(FI)) %*% score
}

b0 <- c(10,10)
b1 <- fs(Y, X[,2], b0)
b2 <- fs(Y, X[,2], b1)
b3 <- fs(Y, X[,2], b2)
b4 <- fs(Y, X[,2], b3)
b5 <- fs(Y, X[,2], b4)

}

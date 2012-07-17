#' .qplot.data
#' This code is taken from qvalue::qplot (from R2.9.0) in order to export the required data for my
#' LimmaGP GenePattern module
#'
#' @param qobj a qvalue object.
#' @param rng see \code{\link[qvalue]{qplot}}
#' @param smooth.df see \code{\link[qvalue]{qplot}}
#' @param smooth.log.pi0 see \code{\link[qvalue]{qplot}}
#' @param \dots unused
#' 
#' @return 
#'	a named list of 4 vectors:
#'	\item{\code{lambda}}{the lambda values (usually 0...0.9)}
#'	\item{\code{pi0}}{the estimated pi0 for each lambda}
#'	\item{\code{spline}}{the cubed spline values for each lambda}
#'	\item{\code{pi0.estimate}}{the final pi0 estimate (from the qobj itself). length 1}
#'  
#' for example:\samp{
#' List of 4
#'  $ lambda      : num [1:20] 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 ...
#'  $ pi0         : num [1:20] 1 0.975 0.961 0.942 0.92 ...
#'  $ spline      : num [1:20] 0.987 0.97 0.953 0.936 0.92 ...
#'  $ pi0.estimate: num 0.793
#' }
#' @author Mark Cowley, 2009-11-30
#' @seealso \code{\link[qvalue]{qplot}}
#' @export
#' @rdname qplot.data
.qplot.data <- function (qobj, rng = c(0, 0.1), smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
	is(qobj, "qvalue") || stop("qobj should be a qvalue object")
	
	q2 <- qobj$qval[order(qobj$pval)]
	if (min(q2) > rng[2]) {
		rng <- c(min(q2), quantile(q2, 0.1))
	}
	
	p2 <- qobj$pval[order(qobj$pval)]

	lambda <- qobj$lambda
	if (length(lambda) == 1) {
		lambda <- seq(0, max(0.9, lambda), 0.05)
	}
	
	pi0 <- rep(0, length(lambda))
	for (i in 1:length(lambda)) {
		pi0[i] <- mean(p2 > lambda[i])/(1 - lambda[i])
	}
	
	if (smooth.log.pi0) 
		pi0 <- log(pi0)
	
	spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
	
	if (smooth.log.pi0) {
		pi0 <- exp(pi0)
		spi0$y <- exp(spi0$y)
	}
	# plot(lambda, pi0, xlab = expression(lambda), ylab = expression(hat(pi)[0](lambda)), pch = ".")
	# lines(spi0)
	res <- list(lambda=lambda, pi0=pi0, spline=spi0$y, pi0.estimate=qobj$pi0)
	return( res )
}
# CHANGELOG:
# 2012-03-07 - changed \verbatim to \samp
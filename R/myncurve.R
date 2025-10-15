#' Displays a shaded curve
#'
#' @param mu integer that represents the mean of a distribution
#' @param sigma integer that represents the standard deviation of a distribution
#' @param a integer
#' @importFrom graphics curve polygon
#' @importFrom stats dnorm pnorm
#'
#' @returns A shaded area between the curve and x axis from negative infinity to x = a.
#'
#' @export
#'
#' @examples
#' myncurve(mu = 10, sigma = 5, a = 6)
#'
#' @name myncurve

utils::globalVariables("x")
myncurve = function(mu, sigma, a){
  curve(dnorm(x,mean=mu,sd=sigma), xlim = c(mu-3*sigma, mu +
                                              3*sigma))
  list(mu = mu, sigma = sigma)

  xcurve = seq(mu - 3*sigma, a, length = 1000)
  ycurve = dnorm(xcurve, mean = mu, sd = sigma)
  polygon(c(mu - 3*sigma, xcurve, a), c(0, ycurve, 0), col = "red")

  prob = pnorm(a, mean = mu, sd = sigma)

  list(mu = mu, sigma = sigma, probability = prob)
}

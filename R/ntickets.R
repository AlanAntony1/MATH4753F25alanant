#' Overbooking Problem
#'
#' @param N The number of seats on the plane.
#' @param gamma The pain probability, aka the probability that we have more passangers show up than seats on the plane.
#' @param p The probability that a passenger shows up.
#' @importFrom stats pbinom uniroot
#' @importFrom graphics abline
#'
#' @returns The amount of seats that an airline should sell, based on the probability that a passenger shows up.
#' @export
#'
#' @examples
#' ntickets(N = 400, gamma = 0.02, p = 0.95)
#'
ntickets <- function(N, gamma, p){
  #Discrete Distribution
  xx <- seq(N, N+40, by = 1) #

  discrete <- sapply(xx, function(n){
    1- gamma - pbinom(N, n, p)
  })

  dIndex <- which.min(abs(discrete))
  nd <- xx[dIndex]

  #Normal Distribution
  normal <- sapply(xx, function(n){
    mean <- n*p
    sd <- sqrt(n*p*(1-p))
    1- gamma - pnorm(N+0.5, mean, sd)
  })

  cRoot <- uniroot(f = function(n){
    mean <- n*p
    sd <- sqrt(n*p*(1-p))
    1 - gamma - pnorm(N+0.5, mean, sd)
  },
  lower = min(xx), upper = max(xx))

  nc <- cRoot$root

  results <- list(
    nd = nd,
    nc = nc,
    N = N,
    p = p,
    gamma = gamma
  )
  plot(xx, discrete, #Data
       ylab = "Objective", xlab = "n", main = paste("Objective vs n to find optimal tickets sold - Discrete", nd, sep = "\n"), #Labels
       type = "b", cex = 0.5, pch = 16, col = "blue") #Lines and Points
  abline(h = 0,
         lty = 1, col = "red")
  abline(v = nd, lty = 1, col = "red")

  plot(xx, normal, #Data
       ylab = "Objective", xlab = "n", main = paste("Objective vs n to find optimal tickets sold - Continuous", nc, sep = "\n"), #Labels
       type = "l", cex = 0.5, pch = 16, col = "blue") #Lines and Points
  abline(h = 0,
         lty = 1, col = "blue")
  abline(v = nc, lty = 1, col = "blue")
  print(results)
  results
}

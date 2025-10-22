#' Overbooking Problem
#'
#' @param N The number of seats on the plane.
#' @param gamma The pain probability, aka the probability that we have more passangers show up than seats on the plane.
#' @param p The probability that a passenger shows up.
#' @importFrom stats pbinom
#'
#' @returns The amount of seats that an airline should sell, based on the probability that a passenger shows up.
#' @export
#'
#' @examples
ntickets <- function(N, gamma, p){
  #Discrete Distribution
  xx <- seq(N, N+200, by = 1)

  discrete <- sapply(xx, function(n){
    pbinom(N, n, p) - 1 + gamma
  })

  index <- which.min(abs(discrete))
  xx[index]
}

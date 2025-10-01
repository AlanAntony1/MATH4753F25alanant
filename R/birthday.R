#' Birthday Problem Functions
#'
#' @param x Total number of people.
#'
#' @returns Returns the chance that at least two people in x number of people have the same birthday.
#' @export
#'
#' @examples
birthday <- function(x){
  1 - exp(lchoose(365,x) + lfactorial(x) - x*log(365))
}

InitWtErgmTerm.Taper <- function(nw, arglist, response=NULL, ...){
  out <- InitErgmTerm.Taper(nw, arglist, response=response, ...)
  out$name <- "wttaper_term"
  out
}

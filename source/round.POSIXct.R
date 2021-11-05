# courtesy of Rui Barradas
# found at https://stat.ethz.ch/pipermail/r-help/2012-June/315336.html
# briefly modified so that the original timezone is preserved
round.POSIXct <- function(x, units = c("mins", "5 mins", "10 mins", "15 mins", "30 mins", "60 mins")){
  if(is.numeric(units)) units <- as.character(units)
  units <- match.arg(units)
  r <- switch(units,
              "mins" = 60,
              "5 mins" = 60*5,
              "10 mins" = 60*10,
              "15 mins" = 60*15,
              "30 mins" = 60*30,
              "60 mins" = 60*60)
  
  H <- as.integer(format(x, "%H"))
  M <- as.integer(format(x, "%M"))
  S <- as.integer(format(x, "%S"))
  tz<-attr(x, "tzone")
  D <- format(x, "%Y-%m-%d")
  D <- as.POSIXct(D, tz=tz)
  secs <- 3600*H + 60*M + S
  return(as.POSIXct(round(secs/r)*r, origin = D, tz = tz))
}
  
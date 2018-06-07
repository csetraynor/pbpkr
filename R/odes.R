#' ODE
#' 
#' @param y 
#'
#' @param times 
#' @param func 
#' @param parms 
#' @param var 
#' @param range_var 
#' @param ... 
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang !!
odes <- function(y, times, func, parms, var,range_var, ...){
  ode_list <- lapply(range_var, function(i) {
    parms[var] = i
    out = deSolve::ode(y = In_Cond, times = times, func = func, parms = parms)
    data.frame(
    x = out[,"time"],
    y = rep(i, length(out)),
    z = out[,"y5"]

    )
  }
  )
  do.call(rbind, ode_list)
}
odes.plot <- function(object){
  object <- unique(object)
  v <- ggplot(object, aes(x, y, z = z))
  v + geom_raster(aes(fill = z))+
    geom_contour(colour = "white")
}






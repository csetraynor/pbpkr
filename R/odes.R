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
    parameters_Norm[var] = i
    out = deSolve::ode(y = In_Cond, times = times, func = PitaODE, parms = parameters_Norm)
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
  v <- ggplot(object, aes(x, y, z = z)) 
  v +stat_density_2d()
}






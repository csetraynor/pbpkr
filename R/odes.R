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
#' 
#' @export
odefun <- function(object, ...){
  UseMethod("odefun") 
}


#' @export
odefun <- function(comp, michaelis){
  #Write ode
  ncomp <- 1:comp
  fun <-  lapply(seq_along(ncomp), function(i){
    
    paste0(ifelse((-1)^(i+1) > 0 , "", "-"),
           'k', i, '*y', i, 
           if(i %in% michaelis ){
             #write michaelis constant
           }) 
  } )
  fun <- as.data.frame(do.call(rbind, fun))
  names(fun) = "fun"
  
  y_var <-  lapply(seq_along(ncomp), function(i){
    
    paste0("dy", i, "dt")
  } )
  y_var <- as.data.frame(do.call(rbind, y_var) )
  names(y_var) = "y_var"
  out <- cbind(y_var, fun)
  
  class(out) <- c("odefun", "data.frame")
  out
}

#' @export
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



# ode_fun <- function(prior = "", class = "", coef = "", group = "", 
#                       resp = "", dpar = "", nlpar = "", bound = "",
#                       ls = list()) {
#   # helper function to create data.frames containing prior information
#   if (length(ls)) {
#     if (is.null(names(ls))) {
#       stop("Argument 'ls' must be named.")
#     }
#     names <- c(
#       "prior", "class", "coef", "group", 
#       "resp", "dpar", "nlpar", "bound"
#     )
#     if (!all(names(ls) %in% names)) {
#       stop("Names of 'ls' must some of ", collapse_comma(names))
#     }
#     for (v in names(ls)) {
#       assign(v, ls[[v]])
#     }
#   }
#   out <- data.frame(
#     prior, class, coef, group, resp, dpar, nlpar, bound, 
#     stringsAsFactors = FALSE
#   )
#   class(out) <- c("brmsprior", "data.frame")
#   out
# }






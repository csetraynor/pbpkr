brmsprior <- function(prior = "", class = "", coef = "", group = "", 
                      resp = "", dpar = "", nlpar = "", bound = "",
                      ls = list()) {
  # helper function to create data.frames containing prior information
  if (length(ls)) {
    if (is.null(names(ls))) {
      stop("Argument 'ls' must be named.")
    }
    names <- c(
      "prior", "class", "coef", "group", 
      "resp", "dpar", "nlpar", "bound"
    )
    if (!all(names(ls) %in% names)) {
      stop("Names of 'ls' must some of ", collapse_comma(names))
    }
    for (v in names(ls)) {
      assign(v, ls[[v]])
    }
  }
  out <- data.frame(
    prior, class, coef, group, resp, dpar, nlpar, bound, 
    stringsAsFactors = FALSE
  )
  class(out) <- c("brmsprior", "data.frame")
  out
}



empty_brmsprior <- function() {
  # define a brmsprior object with zero rows
  char0 <- character(0)
  brmsprior(
    prior = char0, class = char0, coef = char0, 
    group = char0, resp = char0, dpar = char0,
    nlpar = char0, bound = char0
  )
}
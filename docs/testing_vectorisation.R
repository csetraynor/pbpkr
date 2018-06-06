rangeCL_BiPi = seq(from = 165, to = 175, by = 1)


test <- odes(y = In_Cond, times = times,
             func = PitaODE, parms = parameters_Hi,
             var = "CL_BiPi", range_var = rangeCL_BiPi )

out <- lapply(CL_BiPi, function(x) {
  parameters_Norm["CL_BiPi"] = x
  ode(y = In_Cond, times = times, func = PitaODE, parms = parameters_Norm)
}
)


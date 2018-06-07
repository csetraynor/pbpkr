rangeVmP = seq(from = 8e6, to = 11e6, by = .1e6)


test <- odes(y = In_Cond, times = times,
             func = PitaODE, parms = parameters_Norm,
             var = "VmP", range_var = rangeVmP )
p <- odes.plot(test)
p + ggtitle("S-C Plot")


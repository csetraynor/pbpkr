#'ODE function
PitaODE <- function(t, In_Cond, parameters) 
{with(as.list(c(In_Cond, parameters)),{
  
  # pita lag compartment
  dy1dt <- -ktransP*y1   
  # Pita gut compartment (y2)
  dy2dt <- - kaP*y2  + BiTrans*(CL_BiPi*y4/VGaBl) + ktransP*y1   
  # Pita liver extracellular space (y3)
  dy3dt <- (kaP*y2                                                                     
            - Cfpp*y3*(VmP/(KmP+y3) + PdPi) 
            - Qh*(y3 - y5) + y4*PdePi
  )/Vext        
  # Pita liver (y4)
  dy4dt <- (Cfpp*y3*(VmP/(KmP+y3) + PdPi)                                          
            - (CL_BiPi + CL_MePi + PdePi)*y4
  )/VH          
  # Pita plasma (y5)
  dy5dt <- (Qh*(y3 - y5) - Cfpp*CL_ur*y5)/VcP                                                            
  
  list(c(dy1dt,dy2dt,dy3dt,dy4dt,dy5dt)) })  }
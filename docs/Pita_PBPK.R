setwd("\\\\GOBO/User46/u/u1551634/Documents/R")

library(deSolve)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(cowplot)

# Model Based on: Takeuchi et al (2014). DMD 45:726-734, with added biliary excretion into the gut and delay compartment for Tlag
# Clinical data: Prueksaritanont et al. (2014). Br.J.Clin.Pharmacol 78(3): 587-598
Prueksaritanont <- read.csv('Pruesksaritanont_2014_Pita.csv',head=TRUE)

# y1-y4 = pita
# assuming that 100% of the dose is available to be absorbed (Fa=1)

# CfPP = fu.plasma/BL:PL 

####################Function##########################################
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

times <- c(seq(0,600,1)) # min
######################################################################
parameters_Norm <- c(
1000000          ,     # Dose, ng, Prueksaritanont et al. (2014). Br.J.Clin.Pharmacol 78(3): 587-598
ktransP = 0.1    ,     #/min, min gastric emptying time = 10min. Hirano et al (2006) 
kaP = 0.07       ,     # FDA doc, ka=1/(MRT_PO - MRT_IV)=1/(4.55-2.18)=0.42h^-1= 0.007min^-1/ktransp
BiTrans = 0.0618 ,     # /min. Gallbladder emptying rate. Guiastrennec et al (2016). CPT Pmetric 5(12):692-700  
CL_BiPi = 169.63 ,     # Biliary CL, Total CL = 394.5mL/min (FDA)*0.43(fraction in faeces) = 169.63ml/min
#fTP = 0.026      ,    # fu, tissue from RED expt #2
Qh = 1450        ,     # Hepatic blood flow (ml/min)
VGaBl = 30       ,     # Gallbladder volume (mL), Guiastrennec et al (2016)
VmP = 10755828   ,     # ng/min/liver. Model C. 0.171*421.46*139*1800*0.6(fraction hepatocytes) = 10755828
Cfpp = 0.009     ,     # fu.pl/BL:PL = 0.004 (Aus TGA, 2013) / 0.425 (Izumi et al (2017). J.Pharm.Sci 106: 2678-2687 )
KmP = 139.17     ,     # ng.  Model C (12.7*421.46*ftp) = 5352*0.026 = 139.17
PdPi = 683       ,     # P.diff into heps. ml/min. Model C, Pdiff =0.68 of trans up, therefore 2009(scaled*0.68)=1366ml/min /2 
PdePi = 464      ,     # P.diff out of heps. ml/min. Model C 0.34 of Pdif, therefore 1366*0.34 = 464  
Vext = 469       ,     # ml. Watanabe et al (2009). JPET 328:652-662
CL_MePi = 181.44 ,     # ml/min. Fujino et al (2003) 2.52ul/min/mg prot * 40(MPPGL) * 1800 = 181.44
VH = 1690        ,     # ml. vol of liver. Davies and Morris (1993)
CL_ur = 11.84    ,     # ml/min.Urine. FDA = 3% of parent in urine
VcP = 16890            # Vz for pita. ml.Vz/F=226*0.51*1000 (Prueksaritanont et al (2014). ) 
)
In_Cond <- c(y1=parameters_Norm[1],y2=0,y3=0,y4=0,y5=0)
#In_Cond <- c(y2=parameters_Norm[1],y3=0,y4=0,y5=0)


#Vz for pita. FDA = 133.2L. Gives very low plasma values. Total Body Water = 42L 
#HPGL:139, Sohlenius-Sternbeck (2006). ToxInVit 20: 1582-1583. wt of liver,1800: Davies and Morris (1993). Pharm Res 10(7): 1093-1095 
#MPPGL: Zhang et al (2015). Sci Rep 5:17671

out <- ode(y = In_Cond, times = times, func = PitaODE, parms = parameters_Norm)
out<-signif(out,3)
LOQ <- rep(1,length(times))
Sim <- data.frame(out)
Sim$time <- Sim$time/60
Sim[,7]<-LOQ
Sim2 <- data.frame(Sim[-1,])
Sim2$y5 <- Sim2$y5*0.425  #to get from Blood to plasma [], Bl:Pl = 0.425

fig1 <- ggplot(data=Sim2)+
  #geom_line(aes(x=time,y=y4),color="red")+    #y4=liver
  geom_line(aes(x=time,y=y5),color="blue")+   #y5=plasma
  geom_line(aes(x=time,y=V7),linetype=2)+
  scale_y_log10()+
  labs(x="Time (h)",y="[pitavastatin] (ng/ml)")
fig1

fig2 <- fig1 + geom_point(data=Prueksaritanont,aes(x=Time,y=DV))+
  geom_errorbar(data=Prueksaritanont,aes(x=Time,ymin=Neg_err,ymax=Pos_err))
fig2

#############################################################################################
parameters_Low <- c(
  1000000          ,     # Dose, ng, Prueksaritanont et al. (2014). Br.J.Clin.Pharmacol 78(3): 587-598
  ktransP = 0.1    ,     #/min, min gastric emptying time = 10min. Hirano et al (2006) 
  kaP = 0.07,      # FDA doc, ka=1/(MRT_PO - MRT_IV)=1/(4.55-2.18)=0.42h^-1= 0.007min^-1/ktransp
  BiTrans = 0.0618 ,     # /min. Gallbladder emptying rate. Guiastrennec et al (2016). CPT Pmetric 5(12):692-700  
  CL_BiPi = 169.63 ,     # Biliary CL, Total CL = 394.5mL/min (FDA)*0.43(fraction in faeces) = 169.63ml/min
  #fTP = 0.026      ,     # fu, tissue from RED expt #2
  Qh = 1450        ,     # Hepatic blood flow (ml/min)
  VGaBl = 30       ,      # Gallbladder volume (mL), Guiastrennec et al (2016)
  VmP = 8604662   ,     # ng/min/liver. Model C. 0.171*421.46*139*1800*0.6(fraction hepatocytes) = 10755828
  Cfpp = 0.009     ,      # fu.pl/BL:PL = 0.004 (Aus TGA, 2013) / 0.425 (Izumi et al (2017). J.Pharm.Sci 106: 2678-2687 )
  KmP = 172.57 ,       # ng.  Model C (12.7*421.46*ftp) = 5352*0.026 = 139.17
  PdPi = 683,#1366      ,     # P.diff into heps. ml/min. Model C, Pdiff =0.68 of trans up, therefore 2009(scaled*0.68)=1366ml/min /2 
  PdePi = 464      ,     # P.diff out of heps. ml/min. Model C 0.34 of Pdif, therefore 1366*0.34 = 464  
  Vext = 469       ,      # ml. Watanabe et al (2009). JPET 328:652-662
  CL_MePi = 181.44 ,      # ml/min. Fujino et al (2003) 2.52ul/min/mg prot * 40(MPPGL) * 1800 = 181.44
  VH = 1690        ,      # ml. vol of liver. Davies and Morris (1993)
  CL_ur = 11.84    ,     # ml/min.Urine. FDA = 3% of parent in urine
  VcP = 16890            # Vz for pita. ml.Vz/F=226*0.51*1000 (Prueksaritanont et al (2014). ) 
)
out <- ode(y = In_Cond, times = times, func = PitaODE, parms = parameters_Low)
out<-signif(out,3)
Sim <- data.frame(out)
Sim$time <- Sim$time/60
#Sim[,7]<-LOQ
Sim3 <- data.frame(Sim[-1,])*0.425  #to get from Blood to plasma [], Bl:Pl = 0.425

parameters_Hi <- c(
  1000000          ,     # Dose, ng, Prueksaritanont et al. (2014). Br.J.Clin.Pharmacol 78(3): 587-598
  ktransP = 0.1    ,     #/min, min gastric emptying time = 10min. Hirano et al (2006) 
  kaP = 0.07,      # FDA doc, ka=1/(MRT_PO - MRT_IV)=1/(4.55-2.18)=0.42h^-1= 0.007min^-1/ktransp
  BiTrans = 0.0618 ,     # /min. Gallbladder emptying rate. Guiastrennec et al (2016). CPT Pmetric 5(12):692-700  
  CL_BiPi = 169.63 ,     # Biliary CL, Total CL = 394.5mL/min (FDA)*0.43(fraction in faeces) = 169.63ml/min
  #fTP = 0.026      ,     # fu, tissue from RED expt #2
  Qh = 1450        ,     # Hepatic blood flow (ml/min)
  VGaBl = 30       ,      # Gallbladder volume (mL), Guiastrennec et al (2016)
  VmP = 12906993   ,     # ng/min/liver. Model C. 0.171*421.46*139*1800*0.6(fraction hepatocytes) = 10755828
  Cfpp = 0.009     ,      # fu.pl/BL:PL = 0.004 (Aus TGA, 2013) / 0.425 (Izumi et al (2017). J.Pharm.Sci 106: 2678-2687 )
  KmP = 105.77    ,       # ng.  Model C (12.7*421.46*ftp) = 5352*0.026 = 139.17
  PdPi = 683,#1366      ,     # P.diff into heps. ml/min. Model C, Pdiff =0.68 of trans up, therefore 2009(scaled*0.68)=1366ml/min /2 
  PdePi = 464      ,     # P.diff out of heps. ml/min. Model C 0.34 of Pdif, therefore 1366*0.34 = 464  
  Vext = 469       ,      # ml. Watanabe et al (2009). JPET 328:652-662
  CL_MePi = 181.44 ,      # ml/min. Fujino et al (2003) 2.52ul/min/mg prot * 40(MPPGL) * 1800 = 181.44
  VH = 1690        ,      # ml. vol of liver. Davies and Morris (1993)
  CL_ur = 11.84    ,     # ml/min.Urine. FDA = 3% of parent in urine
  VcP = 16890            # Vz for pita. ml.Vz/F=226*0.51*1000 (Prueksaritanont et al (2014). ) 
)

out <- ode(y = In_Cond, times = times, func = PitaODE, parms = parameters_Hi)
out<-signif(out,3)
Sim <- data.frame(out)
Sim$time <- Sim$time/60
#Sim[,7]<-LOQ
Sim4 <- data.frame(Sim[-1,])*0.425  #to get from Blood to plasma [], Bl:Pl = 0.425

theme_set(theme_cowplot(font_size=8,font_family = "Times New Roman"))

fig3 <- fig2 + #geom_line(data=Sim3,aes(x=time,y=y4),color="red",linetype=2)+    #Low CLup
                geom_line(data=Sim3,aes(x=time,y=y5),color="blue",linetype=2)
fig4 <- fig3 + #geom_line(data=Sim4,aes(x=time,y=y4),color="red",linetype=3)+    #Hi CLup
              geom_line(data=Sim4,aes(x=time,y=y5),color="blue",linetype=3)
              
fig4

#ggsave("Pita_PBPK_comp.png",fig4,width=7,height=3.8,units="in",dpi=600) 

# fig3<-ggplot(data=Sim2)+
#   geom_line(aes(x=time,y=y2))+
#   scale_y_log10()
# fig3
  
###########################################################


##################################################################################################################################################

####################Function##########################################
# PitaODE <- function(t, In_Cond, parameters) 
# {with(as.list(c(In_Cond, parameters)),{
#   dy1dt <- -ktrans*y1                                                                     # Pita Lag compartment (y1)
#   dy2dt <- - Fa*kaP*y2 + Fa*CL_BiPi*fTP*y4 + ktrans*y1                                      # Pita Gut (y2)
#   dy3dt <- (Fa*kaP*y2                                                                   # Pita liver Extracellular Space (y3)
#             - Cfpp*((VmP*y3/(KmP+y3)) + y3*PdPi) 
#             - Qh*y3 + Qh*y5 
#   )/Vext        
#   dy4dt <- (Cfpp*(((VmP*y3/(KmP+y3)) + y3*PdPi))                                       # Pita liver (y4)
#             - (CL_BiPi + CL_MePi)*fTP*y4
#   )/VH          
#   dy5dt <- (Qh*y3 - Qh*y5)/VcP                                                          # Pita plasma (y5)
#   
#   list(c(dy1dt,dy2dt,dy3dt,dy4dt,dy5dt)) })  }
# 
# times <- c(seq(0,600,1)) # min


# PiEtODE <- function(t, In_Cond, parameters) 
# 		{with(as.list(c(In_Cond, parameters)),{
# 			# odes for media and cell
# 			dy1dt <- Pita*Fa -kaP*y1 + CL_BiPi*fTP*y3/(1+(fTE*y7/KmBE)) # Pita Gut (y1)
# 			dy2dt <- (kaP*y1 + Qh*y4 - Qh*y2 
# 			          - (VmP*y2*CfPP/(Km*(1+(y6*CfPE)/KiE)+y2*Cfpp)) - y2*CfPP*kfP)/Vext #Pita liver ESpace (y2)
# 			dy3dt <- ((VmP*y2*CfPP/(Km*(1+(y6*CfPT)/KiE)+y2*Cfpp)) + y2*CfPP*PdifPi 
# 			          - CL_BiPi*fTP*y3/(1+(fTE*y7/KmBE)) - CL_MePi*fTP*y3)/VH # Pita liver (y3)
# 			dy4dt <- (Qh*y2 - Qh*y4)/VcP
# 			
# 			# return ode solutions
# 			list(c(dy1dt,dy2dt,dy3dt,dy4dt)) })  }

#times <- seq(0, 30, by = 0.25)  
  


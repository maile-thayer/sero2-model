###############################################################################
# Run the 4-serotype SEIR-SEI dengue model.
# ~~~
# 2024-06-03
###############################################################################

library('deSolve')


# Load necessary functions.
source('Garcia-Carreras model/seirsei.R')


# Time points at which we want estimates (in days).
times = seq(from = 1, to = 50 * 365, by = 7)

# Duration of cross-protection between serotypes (in months).
cross_prot = 12

# Vector species (choice between "Ae. aegypti" and "Ae. albopictus").
species = 'Ae. aegypti'



###############################################################################
# Temperature functions to be used in the model.
###############################################################################

# This can be any function that returns a temperature when given a time point
# (in days). I provide two examples here:
# - Sinusoidal temperature (idealisation of annual seasonality).
# - Function that interpolates a provided (in this case made up) weekly dataset.


#' Sinusoidal temperature function.
#'
#' @param tp :numeric: time point (in days).
#' @param t_amp :numeric: amplitude (mean to peak).
#' @param t_mean :numeric: mean temperature.
#' @param t_seas :numeric: day on which the function peaks.
#'
#' @return an interpolating function.

TempFunSine = function(tp, t_amp = 5, t_mean = 28, t_seas = 0)
{
  t_s = (t_seas / 365.25 - 0.25) * 2 * pi
  return(t_amp / 2 * sin(2 * pi * tp / 365.25 - t_s) + t_mean)
}


# Made up weekly temperature data.
x = 1:(52 * 5)
dat = data.frame(times = x, 
                 temps = 10 / 2 * sin(2 * pi * x / 52) + 26 + 
                   runif(n = length(x), min = -5, max = 5))

TempFunDat = approxfun(dat$temps, x = dat$times * 7)



################################################################################
# Assemble parameters to be used in the model (in units of /day where
# pertinent).
################################################################################

pars = list(N_h = 1e6,          # Number of hosts.
            M = 1.5,            # Ratio of vectors to hosts.
            species = species,  # Vector species.
            r_h = 3.55e-5,      # Host birth rate.
            mu_h = 3.55e-5,     # Host mortality rate.
            E_a = 0.05,         # Activation energy of carrying capacity (see note in seirsei.R).
            delta_h = 1/5.9,    # Incubation rate.
            eta_h = 1/5,        # Recovery rate.
            xi_h = 1 / (cross_prot / 12 * 365.25),  # Inverse of cross-protection.
            TempFun = TempFunSine,  # Temperature function to be used.
            times = times,      # Time points at which to produce estimates.
            trickle = 1e-5)     # Added to vector compartments (see note in seirsei.R).




# Starting conditions.
# ~~~
# This generates a vector with the starting numbers in each of the compartments.
# Note that parameter `start_conds` is a choice between two options: 
# - 'susc' which assumes all hosts but 10 are susceptible, the remaining 10
#   allocated to the first-infection compartments of each serotype (1, 2, 3, and
#   4 individuals randomly allocated to each serotype). 
# - 'random', in which all individuals are randomly allocated across all
#   compartments.
xstart = GetStartConds(pop_tots = pars$N_h, M = pars$M, start_conds = 'susc')

dxdt = vector('numeric', length = length(xstart))
names(dxdt) = paste0('d', names(xstart))

z = ode(func = RunSEIRSEI4Sero, times = pars$times, y = xstart, parms = pars,
        rtol = 1e-12)

require(tidyverse)
z2 <- as.data.frame(z) %>%
  mutate(Eh_all = Eh1+Eh2+Eh3+Eh4+
            Eh12+Eh21+Eh31+Eh41+
            Eh13+Eh23+Eh32+Eh42+
            Eh14+Eh24+Eh34+Eh43,
         Ih_all = Ih1+Ih2+Ih3+Ih4+
           Ih12+Ih21+Ih31+Ih41+
           Ih13+Ih23+Ih32+Ih42+
           Ih14+Ih24+Ih34+Ih43,
         Ph_all = Ph1+Ph2+Ph3+Ph4,
         Rh_all = Rh1+Rh2+Rh3+Rh4)
plot(z2$time,z2$Eh_all, type="l",ylim=c(0,1500))

newI <- mutate(as_tibble(z),
  inf1 = Ih1 + Ih21 + Ih31 + Ih41,
  inf2 = Ih2 + Ih12 + Ih32 + Ih42,
  inf3 = Ih3 + Ih13 + Ih23 + Ih43,
  inf4 = Ih4 + Ih14 + Ih24 + Ih34,
  inf_all = inf1 + inf2 + inf3 + inf4) %>%
  select(time, inf1, inf2, inf3, inf4, inf_all) %>%
  slice(-(1:(52*3)))
plot(1, 1, type='n', xlim=range(newI$time), ylim=range(select(newI, -time)))
lines(select(newI, time, inf_all))
lines(select(newI, time, inf1), col='blue')
lines(select(newI, time, inf2), col='red')
lines(select(newI, time, inf3), col='darkgreen')
lines(select(newI, time, inf4), col='darkorange')

##### STACKED AREA PLOTS FOR COMPARTMENTS
source("data_stackedareaplot-MAJ.R")
result <- as.data.frame(z)
data <- data_stackedareaplot2(list(
  Sh_0 = result$Sh, 
  Sh_secondary = result$Rh1 + result$Rh2 + result$Rh3 + result$Rh4,
  EI_primary = result$Eh1 + result$Eh2 + result$Eh3 + result$Eh4 + 
    result$Ih1 + result$Ih2 + result$Ih3 + result$Ih4,
  EI_secondary = result$Eh12 + result$Eh13 + result$Eh14 + 
    result$Eh21 + result$Eh23 + result$Eh24 + 
    result$Eh31 + result$Eh32 + result$Eh34 +  
    result$Eh41 + result$Eh42 + result$Eh43 +  
    result$Ih12 + result$Ih13 + result$Ih14 + 
    result$Ih21 + result$Ih23 + result$Ih24 + 
    result$Ih31 + result$Ih32 + result$Ih34 +  
    result$Ih41 + result$Ih42 + result$Ih43,
  R_primary = result$Ph1 + result$Ph2 + result$Ph3 + result$Ph4,
  R_secondary = result$Rh
))


S.col <- brewer.pal(9, "GnBu")[8:9]
E.col <- brewer.pal(9, "Purples")[8:9]
I.col <- brewer.pal(9, "YlOrRd")[8:9]
R.col <- brewer.pal(9, "BuGn")[8:9]

# Population distribution of each compartment over timeframe
ggplot(data, aes(x = time, y = value_perc, fill = group)) +
  geom_area() +
  theme_classic() +
  ggtitle("4-serotype model") +
  scale_fill_manual(values = c(S.col, I.col, R.col))


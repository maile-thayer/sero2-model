###############################################################################
# Run the 4-serotype SEIR-SEI dengue model.
# ~~~
# 2024-06-03
###############################################################################

library('deSolve')


# Load necessary functions.
source('seirsei.R')


# Time points at which we want estimates (in days).
times = seq(from = 1, to = 5000, by = 7)

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


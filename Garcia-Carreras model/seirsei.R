###############################################################################
# 4-serotype SEIR-SEI dengue model definitions and equations.
# ~~~
# 2024-06-03
###############################################################################

library('deSolve')



###############################################################################
# Functions for temperature-dependent rates.
# ~~~
# These are rates given in the SI for Mordecai et al. (2017; doi:
# 10.1371/journal.pntd.0005568). The parameter estimates for both Ae. aegypti
# and Ae. albopictus are given, for convenience, in temperature_pars.csv, but
# the parameter values are then hardcoded in the functions below for the sake of
# efficiency. The SI comes with 95% CI intervals, but not with a full posterior
# (I assume there's one), which would be needed to leverage the uncertainty if
# necessary. 
###############################################################################

# General forms of the temperature functions.

# Briere function.
Briere = function(x, a, Tmin, Tmax)
{
  if ((x >= Tmin) & (x <= Tmax) & is.finite(x))
    return(a * x * (x - Tmin) * sqrt(Tmax - x))
  else
    return(0)
}


# Quadratic function.
Quad = function(x, a, Tmin, Tmax)
{
  if((x >= Tmin) & (x <= Tmax) & is.finite(x))
    return(a * (x - Tmin) * (x - Tmax))
  else
    return(0)
}



# Temperature-dependent rates.

# Eggs laid per femail per day.
EggsLaid = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = 8.56e-3
    Tmin = 14.58
    Tmax = 34.61
  } else if (sp == 'Ae. albopictus')
  {
    a = 4.88e-2
    Tmin = 8.02
    Tmax = 35.65
  }
  
  return(Briere(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Mosquito egg-to-adult survival probability.
SurvProb = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = -5.99e-3
    Tmin = 13.56
    Tmax = 38.29
  } else if (sp == 'Ae. albopictus')
  {
    a = -3.61e-3
    Tmin = 9.04
    Tmax = 39.33
  }

  return(Quad(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Mosquito egg-to-adult development rate per day.
DevlRate = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = 7.86e-5
    Tmin = 11.36
    Tmax = 39.17
  } else if (sp == 'Ae. albopictus')
  {
    a = 6.38e-5
    Tmin = 8.60
    Tmax = 39.66
  }

  return(Briere(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Mosquito mortality rate per day (inverse of lifespan).
MortRate = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = -1.48e-1
    Tmin = 9.16
    Tmax = 37.73
  } else if (sp == 'Ae. albopictus')
  {
    a = -1.43
    Tmin = 13.41
    Tmax = 31.51
  }

  return(1 / Quad(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Biting rate per day.
BiteRate = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = 2.02e-4
    Tmin = 13.35
    Tmax = 40.08
  } else if (sp == 'Ae. albopictus')
  {
    a = 1.93e-4
    Tmin = 10.25
    Tmax = 38.32
  }

  return(Briere(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Probability of mosquito infection.
ProbInfc = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = 4.91e-4
    Tmin = 12.22
    Tmax = 37.46
  } else if (sp == 'Ae. albopictus')
  {
    a = 4.39e-4
    Tmin = 3.62
    Tmax = 36.82
  }

  return(Briere(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Probability of mosquito infectiousness.
ProbInfs = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = 8.49e-4
    Tmin = 17.05
    Tmax = 35.83
  } else if (sp == 'Ae. albopictus')
  {
    a = 7.35e-4
    Tmin = 15.84
    Tmax = 36.40
  }

  return(Briere(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}


# Virus extrinsic incubation rate per day.
IncuRate = function(tmp, sp)
{
  if (sp == 'Ae. aegypti')
  {
    a = 6.65e-5
    Tmin = 10.68
    Tmax = 45.90
  } else if (sp == 'Ae. albopictus')
  {
    a = 1.09e-4
    Tmin = 10.39
    Tmax = 43.05
  }

  return(Briere(x = tmp, a = a, Tmin = Tmin, Tmax = Tmax))
}



#' Carrying capacity.
#'
#' @param x :numeric: temperature (in Celsius).
#' @param Tref :numeric: reference temperature (temperature at which carrying
#'        capacity is maximised; in Celsius).
#' @param N :numeric: maximum carrying capacity (number of mosquitoes).
#' @param Ea :numeric: activation energy of the carrying capacity temperature
#'        function (the sensitivity of carrying capacity to temperature).
#' @param sp :character: mosquito species (choice between "Ae. aegypti" and 
#'        "Ae. albopictus").
#'
#' @return numeric with carrying capacity for the given parameters.
#'
#' NOTE: The carrying capacity term is not particularly well justified, but is
#' necessary to cap the vector population. See SI section on model caveats in
#' Garcia-Carreras et al. (2022; doi: 0.1371/journal.pbio.3001160) for further
#' details.
#'
#' NOTE: The activation energy term `Ea`, which is the degree of temperature
#' sensitivity of the carrying capacity, is ambiguous in Huber et al. (2018);
#' the text mentions one value, while the SI code mentions another (0.5 and
#' 0.05, respectively). It would perhaps be worth keeping this in mind when
#' performing robustness analyses.
CarrCapFun = function(x, Tref, N, E_a, sp)
{
  k = 8.617e-5  # Boltzmann-Arrhenius constant.

  term = EggsLaid(Tref, sp) * SurvProb(Tref, sp) * DevlRate(Tref, sp) / MortRate(Tref, sp)

  return((term - MortRate(Tref, sp)) / term *
         N * exp((-E_a * (x - Tref)^2) / (k * (x + 273.15) * (Tref + 273.15))))
}

CarrCap = function(tmp, N_max, sp)
{
  return(CarrCapFun(x = tmp, Tref = 29, N = N_max, E_a = 0.05, sp = sp))
}



###############################################################################
# Starting conditions.
###############################################################################

#' Generate starting conditions for the model.
#'
#' @param pop_tots :numeric: Total human population size.
#' @param M :numeric: Ratio of mosquitoes to humans.
#' @param start_conds :character: Choice between 'random' and 'susc' - where to
#'        allocate individuals. A choice of:
#'        - 'random': pop_tots are allocated across all compartments at random;
#'        - 'susc': all hosts are assumed to be susceptible, except for 1, 2, 3,
#'        and 4 individuals that are assigned one each to the infected serotype
#'        compartments.

GetStartConds = function(pop_tots, M, start_conds = 'random')
{
  # Hosts.
  hs  = expand.grid(i1 = 1:4, i2 = 1:4)
  hs  = hs[-which(hs$i1 == hs$i2), ]
  hs1 = paste0(rep(c('Eh', 'Ih', 'Ph', 'Rh'), each = 4), rep(1:4, 4))
  hs2 = paste0(hs$i1, hs$i2)
  hs2 = paste0(rep(c('Eh', 'Ih'), each = length(hs2)), rep(hs2, 2))
  hs  = c('Sh', hs1, hs2, 'Rh')

  # Vectors.
  vs = c('Sv', paste0('Ev', 1:4), paste0('Iv', 1:4))

  comps = c(hs, vs)

  xstart = rep(0, length(comps))

  if (start_conds == 'random')
  {
    GetRandStart = function(x)
    {
      rand_start = runif(length(which(x))) 
      return(rand_start / sum(rand_start))
    }

    g = grepl('h', comps)
    xstart[g] = GetRandStart(g) * pop_tots
    g = grepl('v', comps)
    xstart[g] = GetRandStart(g) * pop_tots * M
  } else if (start_conds == 'susc')
  {
    ghs = grep('Sh', comps)
    ghi = grep('Ih[1-4]$', comps)
    gvs = grep('Sv', comps)
    gvi = grep('Iv[1-4]', comps)

    xstart[ghs] = pop_tots - 10
    xstart[ghi] = sample.int(4)
    xstart[gvs] = pop_tots * M 
  } 

  names(xstart) = comps

  return(xstart)
}




###############################################################################
# SEIR-SEI 4-serotype dengue model.
# ~~~
# This model is based on that by Huber et al. (2018; doi:
# 10.1371/journal.pntd.0006451), but extended to four serotypes (rather than the
# one in their paper), including complete immunity after two infections, and
# cross-protection between serotypes, but with no antibody-dependent
# enhancement. For further details see Garcia-Carreras et al (2022; doi:
# 10.1371/journal.pbio.3001160), and in particular, the section on caveats in
# the corresponding SI.
###############################################################################

RunSEIRSEI4Sero = function(tms, x, pars)
{
  with(as.list(c(pars)), {
         seros      = seq_len(4)  # 4 serotypes.
         comp_names = names(x)    # Get the names of the compartments.

         # Temperature at current time.
         tmp = TempFun(tms)

         # Temperature dependent rates (example rates at 26 degrees). 
         biterate = BiteRate(tmp, sp = species)  # = 0.2492967
         probinfs = ProbInfs(tmp, sp = species)  # = 0.6194137
         probinfc = ProbInfc(tmp, sp = species)  # = 0.59552
         mortrate = MortRate(tmp, sp = species)  # = 0.034206
         incurate = IncuRate(tmp, sp = species)  # = 0.1181627
         eggslaid = EggsLaid(tmp, sp = species)  # = 7.45787
         survprob = SurvProb(tmp, sp = species)  # = 0.9157967
         devlrate = DevlRate(tmp, sp = species)  # = 0.108575
         carrcapa = CarrCap(tmp = tmp, N_max = N_h * M, sp = species)  # = 188464.8

         # Keep track of which compartments are vectors.
         vs  = grep('v', comp_names)
         N_v = sum(x[vs])



         #######################################################################
         # Equations for hosts.
         #######################################################################

         # Hosts infected by vectors, per serotype.
         hivi = biterate * probinfs * x[paste0('Iv', seros)] / N_h

         # Hosts susceptible.
         dxdt['dSh'] = r_h * N_h -  # Birth rate.
           (sum(hivi) + mu_h) * x['Sh']  # Infection + death.


         for (se in seros)
         {
           # Hosts exposed.
           dxdt[paste0('dEh', se)] = hivi[se] * x[paste0('Sh')] - 
             (delta_h + mu_h) * x[paste0('Eh', se)]
         }

         for (se in seros)
         {
           # Hosts infected.
           dxdt[paste0('dIh', se)] = delta_h * x[paste0('Eh', se)] - 
             (eta_h + mu_h) * x[paste0('Ih', se)]
         }

         for (se in seros)
         {
           # Hosts cross-protected.
           dxdt[paste0('dPh', se)] = eta_h * x[paste0('Ih', se)] -
             (xi_h + mu_h) * x[paste0('Ph', se)]
         }

         for (se in seros)
         {
           # Hosts recovered from first infection.
           dxdt[paste0('dRh', se)] = xi_h * x[paste0('Ph', se)] - 
             (sum(hivi[-se]) + mu_h) * x[paste0('Rh', se)] 
         }

         for (se1 in seros)
         {
           for (se2 in seros)
           {
             if (se1 != se2)
             {
               # Hosts exposed the second time.
               dxdt[paste0('dEh', se1, se2)] = hivi[se2] *
                 x[paste0('Rh', se1)] - 
                 (delta_h + mu_h) * x[paste0('Eh', se1, se2)]
             }
           }
         }

         for (se1 in seros)
         {
           for (se2 in seros)
           {
             if (se1 != se2)
             {
               # Hosts with second infection.
               dxdt[paste0('dIh', se1, se2)] = delta_h * 
                 x[paste0('Eh', se1, se2)] -
                 (eta_h + mu_h) * x[paste0('Ih', se1, se2)] 
             }
           }
         }

         # Hosts fully immune.
         dxdt['dRh'] = eta_h * 
           sum(x[grep('Ih[1-4][1-4]', comp_names)]) -
           mu_h * x['Rh']



         #######################################################################
         # Equations for vectors.
         #######################################################################

         # Mosquito birth rate.
         birth_rate = eggslaid * survprob * devlrate * 
           N_v * (1 - N_v / carrcapa) / mortrate 


         vihi = vector('numeric', length = 4)

         for (se in seros)
         {
           # Vectors infected by hosts.
           vihi[se] = biterate * probinfc * 
             (x[paste0('Ih', se)] + sum(x[paste0('Ih', seros[-se], se)])) / N_h
         }


         # Vector susceptibles.
         dxdt['dSv'] = birth_rate - (sum(vihi) + mortrate) * x['Sv'] 

         for (se in seros)
         {
           # Vectors exposed.
           dxdt[paste0('dEv', se)] = vihi[se] * x['Sv'] - 
             (incurate + mortrate) * x[paste0('Ev', se)]
         }

         for (se in seros)
         {
           # Vector infecteds.
           #
           # NOTE: the `trickle` term is there to prevent hosts and vectors
           # becoming very very low and taking a long time to recover.
           dxdt[paste0('dIv', se)] = incurate * 
             x[paste0('Ev', se)] -
             mortrate * x[paste0('Iv', se)] + trickle
         }

         return(list(dxdt, beta_h = biterate * probinfs, temp = tmp))
       })
}


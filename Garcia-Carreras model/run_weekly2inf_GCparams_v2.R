rm(list=ls())
###############################################################
###############################################################
##### 1. Install JuliaCall in R, load packages and environment
###############################################################
###############################################################

# DOES NOT WORK ON R VERSION 4.4.0+

# install.packages("JuliaCall")
library(ggplot2)
library(scales)
library(RColorBrewer)
library(JuliaCall)
library(LaplacesDemon)

# load temperature functions from Bernardo
source('Garcia-Carreras model/seirsei.R')


# install_julia() #DONT RUN UNLESS YOU DON'T HAVE JULIA INSTALLED
#julia <- julia_setup(JULIA_HOME = "C:/Users/ruu6/AppData/Local/Programs/Julia-1.7.2/bin")
#julia <- julia_setup(JULIA_HOME = '~/.julia/juliaup/julia-1.10.5+0.aarch64.apple.darwin14/bin')
julia <- julia_setup(JULIA_HOME = '/Users/michael/.julia/juliaup/julia-1.10.5+0.aarch64.apple.darwin14/bin')

############################################
# setwd("~/GitHub/sero2-model/") #set working directory to github project repo if you haven't already
# julia_install_package("Random")
# julia_install_package_if_needed("Random")
julia_library("Random, Distributions, DataFrames, LinearAlgebra, CSV, JLD2, Dates")

#julia_source("1_dengue_4st_imp.jl") # load main model function
#julia_source("2_run_4st_model_sims_randinit_imp.jl") # load function to run main model

julia_source("1_dengue_2inf.jl")
#x <- julia_eval("dengue_2inf!(u0[1], par, 2, 7)")
julia_source("2_run_2inf_model.jl")

julia_exists("dengue_2inf!") # test whether functions loaded and exist :)
julia_exists("multinom_samp")
julia_exists("run_2inf_model!")

##############################################################
###############################################################
##### 2. Initialize parameters and populations/compartments
###############################################################
###############################################################

##### PARAMETERS
# set parameters using Wearing et al paper https://www.pnas.org/doi/full/10.1073/pnas.0602960103#:~:text=Rapid%20and%20unplanned%20urbanization%2C%20increased,contributed%20to%20dramatic%20dengue%20resurgence.


# set timeframe
dt <- julia_eval("dt = 7;") # time step for the model in days (all parameters in days)
tmax <- julia_eval("tmax = 365 * 50;") # in days

# set # simulations to run-- doing it in R
# nsims=2
nsims <- julia_eval("nsims = 3;")

#detach(model_parameters)
######## HUMAN PARAMETERS ########
pop <- 1e6 #3255603

model_parameters <- list(
  # human
  bh = 3.55e-5, #1/(60*365) * 7, # birthrate; from paper
  mu_h = 3.55e-5, #1/(60*365) * 7, # human death rate; set equal to birthrate
  m = 1.5, # initial ratio of mosquitoes to humans; from literature
  c = 0.5, # contact rate (biting rate); from paper
  p_h = 0.3, # probability of a human being infected by an infected mosquito bite; paper 0.38
  #beta_h = rep(0.5 * 0.38, length(tmax)), # seasonal transmission- mosquito to human
  # beta_h = 0.5 * 0.38 * (0.3 * cos((2*pi*(1:(tmax+365)) + 5.295594)/365) + 1), #seasonal transmission- mosquito to human
  
  p_IIP = 1/5.9, # progression rate out of human E state; from paper
  p_IP = 1/5, # progression rate out of human I state; from paper
  p_R = 1/365, # progression rate out of human cross-protective state; currently at 1 year; from paper
  # p_imp = rep((100/365 * 7)/16,16),
  
  #mosquito
  bm = 1/14, # mosquito egg laying rate ("birth rate"); from paper
  phi_m = 1, # emergence rate; set equal to 1 to stabilize model dynamics
  p_m = 0.3, # probability of mosquito being infected by biting infectious human; set equal to p_h like paper
  # beta_m = 0.5 * 0.38, # human-to-mosquito transmission rate
  mu_m = 1/14, # adult mosquito mortality rate, set equal to egg laying rate
  mu_mL = 0.0, # larval mosquito death rate; initially set to 0
  p_EIP = 1/10, # progression rate out of mosquito E state; from paper
  
  daily_import = 1 # total across serotypes
)

attach(model_parameters)

for (i in 1:length(model_parameters)) {
  julia_assign(names(model_parameters)[i], model_parameters[[i]])
}

# Temperature dependent rates (example rates at 26 degrees). 
species <- 'Ae. aegypti'
# Temperature at current time.
TempFunSine = function(tp, t_amp = 5, t_mean = 28, t_seas = 0) {
  t_s = (t_seas / 52 - 0.25) * 2 * pi
  return(t_amp / 2 * sin(2 * pi * tp / 365 - t_s) + t_mean)
}
tmp = TempFunSine(0:(tmax+365))
biterate = Vectorize(BiteRate, 'tmp')(tmp, sp = species)  # = 0.2492967
probinfs = Vectorize(ProbInfs, 'tmp')(tmp, sp = species)  # = 0.6194137
probinfc = Vectorize(ProbInfc, 'tmp')(tmp, sp = species)  # = 0.59552

beta_h <- biterate * probinfs
#julia_assign('beta_h', beta_h)
beta_h_orig <- c * p_h * (0.3 * cos((2*pi*(0:(tmax+365)) + 5.295594) / 365) + 1) #seasonal transmission- mosquito to human
julia_assign('beta_h', 1.8 * beta_h_orig)

beta_m <- biterate * probinfc
#julia_assign('beta_m', beta_m)
beta_m_orig <- c * p_m # human-to-mosquito transmission rate
julia_assign('beta_m', rep(beta_m_orig, length(beta_h)))

(m * median(beta_h) * median(beta_m) * p_EIP) / (mu_m * (p_IP + bh) * (p_EIP + mu_m))
(m * c * p_h * c * p_m * p_EIP) / (mu_m * (p_IP + bh) * (p_EIP + mu_m))

m*c*p_m*1/p_IP*exp(-mu_m/p_EIP) * c*p_h*1/mu_m

c * p_m * 100 #lambdam
c * p_m * 100 * m*pop * 1/p_IP

# save all params into one vector to input into model function
julia_command("par = (; bh, beta_h, bm, phi_m, beta_m, mu_m, mu_mL, mu_h, p_IIP, p_IP, p_EIP, p_R, m, daily_import)") 

######################################################
########### COMPARTMENTS
julia_assign("pop", pop) # Initial total population of PR

compartment_names <- c('Sh_0', 'Eh_prim', 'Ih_prim', 'Rh_prim', 
  'Sh_sec', 'Eh_sec', 'Ih_sec', 'Rh_sec', 
  'Ih_imp',
  'Lm', 'Sm', 'Em', 'Im')
init_compartments <- matrix(0, nrow=length(compartment_names), ncol=nsims, 
  dimnames = list(compartment_names, paste0('V', 1:nsims)))

# initial_s <- pop
# initial_s2 <- rep(0, 4) #rep(0.2/4, 4)
# initial_r1 <- rep(0, 4) #rep(0.2/4, 4)
# initial_r2 <- rep(0, 12) #rep(0.1/12, 12)
# #init_h <- round(pop * t(rdirichlet(nsims, c(initial_s, initial_s2, initial_r1, initial_r2))))
# # rdirichlet too variable, Sh_0 should always be high
# init_h <- rmultinom(nsims, size=pop, prob=c(initial_s, initial_s2, initial_r1, initial_r2))
# rownames(init_h) <- c('Sh_0', 'Sh_sec',  
#   "Rh_prim", "Rh_sec", "Rh_3", "Rh_4", "Rh_12", 'Rh_13', 'Rh_14', 'Rh_21', 'Rh_23', 'Rh_24', 
#   'Rh_31', 'Rh_32', 'Rh_34', 'Rh_41', 'Rh_42', 'Rh_43')

#init_compartments[rownames(init_h), ] <- init_h
init_compartments["Sh_0", ] <- pop

# initialize mosquitos
init_compartments["Sm", ] <- pop * model_parameters$m

# use imported cases to initialize transmission
# for (imp_compartment in c('Ih_imp_1', 'Ih_imp_2', 'Ih_imp_3', 'Ih_imp_4')) {
# #  init_compartments[imp_compartment, ] <- rpois(nsims, 10)
#   init_compartments[imp_compartment, ] <- rnbinom(nsims, mu=10, size = 0.1)
# }

for (i in 1:nsims) {
  #init_compartments[c('Ih_imp_1', 'Ih_imp_2', 'Ih_imp_3', 'Ih_imp_4'), i] <- 
  #  sample.int(4)
  init_compartments['Ih_imp', i] <- rpois(1, 100)
}


# transfer initial values to Julia
julia_assign("init_compartments", init_compartments)
julia_assign("input_compartment_names", compartment_names)
julia_command("input_compartment_names = Tuple(Symbol.(input_compartment_names))")
julia_command("u0 = [];")
for (i in 1:nsims) {
  julia_command(paste0('push!(u0, NamedTuple{input_compartment_names}(init_compartments[:, ', i, ']))'))
}

######
# Vector of all compartments and sims
julia_command("outcomes0 = [];")
julia_command("for i=1:nsims
                push!(outcomes0, (; hpop = pop, mpop = pop * m, lpop = 0, 
                    hbirths = 0, hdeaths= 0, newcases_all_h = 0, newcases_all_m = 0, 
                    p_infect_h = 0, p_infect_m = 0, 
                    newcases_primary_h = 0, newcases_secondary_h = 0, tot_imp = 0))
              end")
            

# list of all initial items of interest AT TIME 0
# 1. vector of individuals in each compartment
# 2. new human infections
# 3. new mosquito infections
# 4. total human population
# 5. total mosquito population
# 6. total human births
# 7. total human deaths
# 8. human force of infection st1
# 9. human force of infection st2
# 10. human force of infection st3
# 11. human force of infection st4
# 12. mosquito force of infection st1
# 13. mosquito force of infection st2
# 14. mosquito force of infection st3
# 15. mosquito force of infection st4
# 16. all new FIRST infections
# X 17. aging (by age) ---- not using in this model
# X 18. deaths by age ---- not using in this model
# 19. importations ---- not using in this model
# 20. new secondary human cases
# 21. total larval population
# 22. new human serotype 1 infections
# 23. new human serotype 2 infections
# 24. new human serotype 2 infections
# 25. new human serotype 2 infections

###############################################################
###############################################################
##### 3. Run model simulations and store output
###############################################################
###############################################################

# set seed for reproducibility
julia_command("Random.seed!(12345);")

result <- julia_eval("run_2inf_model!(nsims, tmax, u0, outcomes0, par, dt)")

###############################################################
###############################################################
##### 4. Plot output
###############################################################
###############################################################
colors <- brewer.pal(8, name = "Dark2")

############### NEWCASES
# all (humans)
plot_all_humans <- function() {
  df <- as.data.frame(result$newcases_all_h)
  df$median <- apply(df, 1, median)
  plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Days", bty = "n")
  for (i in 1:nsims) {
    lines(result$time, df[, i], col = alpha(colors[1], 0.3))
  }
  lines(result$time, df$median, col = "black")

  for (i in 1:nsims) {
    y = df[, i]
    lw1 <- loess(y ~ result$time)
    # j <- order(df[, i])
    lines(result$time, lw1$fitted, col="red")
  }
}
plot_all_humans()

# primary and secondary
plot_humans_by_type <- function(sims_to_plot) {
  df1 <- as.data.frame(result$newcases_primary_h)
  df2 <- as.data.frame(result$newcases_secondary_h)
  
  plot(NA, NA, 
    xlim = c(0, tmax), 
    ylim = c(1, 500), 
#    ylim = c(1, max(df1, df2)), 
    #log='y',
    ylab = "New infections", xlab = "Days", bty = "n")
  for (i in 1:sims_to_plot) {
    lines(result$time, df1[, i], col = alpha(colors[1], 0.3))
    lines(result$time, df2[, i], col = alpha(colors[2], 0.3))
  }
}
plot_humans_by_type(nsims)
plot_humans_by_type(1)


############### FOI
# serotype 1 (humans)
plot_foi <- function() {
  df <- as.data.frame(result['p_infect_h'])
  df$median <- apply(df, 1, median)
  plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), 
    ylab = "FOI", xlab = "Days", bty = "n")
  for (i in 1:nsims) {
    lines(result$time, df[, i], col = alpha(colors[4], 0.15))
  }
  lines(result$time, df$median, col = colors[4])
}
plot_foi()

##### STACKED AREA PLOTS FOR COMPARTMENTS
source("data_stackedareaplot-MAJ.R")
data <- data_stackedareaplot2(list(
    Sh_0 = result$Sh_0, 
    Sh_secondary = result$Sh_sec,
    EI_primary = result$Eh_prim + result$Ih_prim,
    EI_secondary = result$Ih_sec + result$Eh_sec,
    R_primary = result$Rh_prim,
    R_secondary = result$Rh_sec
  ), times = result$time)
  
  
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


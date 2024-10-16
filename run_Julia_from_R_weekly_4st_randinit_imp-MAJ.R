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

julia_source("1_dengue_4st_imp-MAJ.jl")
#x <- julia_eval("dengue_4st_imp!(u0all[1], par, 2)")
julia_source("2_run_4st_model_sims_randinit_imp_MAJ.jl")

julia_exists("dengue_4st_imp!") # test whether functions loaded and exist :)
julia_exists("multinom_samp")
julia_exists("run_4st_model_sims_randinit_imp!")

##############################################################
###############################################################
##### 2. Initialize parameters and populations/compartments
###############################################################
###############################################################

##### PARAMETERS
# set parameters using Wearing et al paper https://www.pnas.org/doi/full/10.1073/pnas.0602960103#:~:text=Rapid%20and%20unplanned%20urbanization%2C%20increased,contributed%20to%20dramatic%20dengue%20resurgence.


# set timeframe
tmax <- julia_eval("tmax = 52*10;")
julia_command("tspan = collect(0.0:(tmax+52));")

# set # simulations to run-- doing it in R
# nsims=2
nsims <- julia_eval("nsims = 3;")


detach(model_parameters)
######## HUMAN PARAMETERS ########
model_parameters <- list(
  # human
  bh = 1/(60*365) * 7, # birthrate; from paper
  m = 2, # initial ratio of mosquitoes to humans; from literature
  c = 0.5 * 7, # contact rate (biting rate); from paper
  p_h = 0.3, # probability of a human being infected by an infected mosquito bite; paper 0.38
  #beta_h = rep(0.5 * 0.38, length(tmax)), # seasonal transmission- mosquito to human
  # beta_h = 0.5 * 0.38 * (0.3 * cos((2*pi*(1:(tmax+365)) + 5.295594)/365) + 1), #seasonal transmission- mosquito to human
  mu_h = 1/(60*365) * 7, # human death rate; set equal to birthrate
  p_IIP = 1/5 * 7, # progression rate out of human E state; from paper
  p_IP = 1/6 * 7, # progression rate out of human I state; from paper
  p_R = 1 * 1/365 * 7, # progression rate out of human cross-protective state; currently at 1 year; from paper
  # p_imp = rep((100/365 * 7)/16,16),
  
  #mosquito
  bm = 1/14 * 7, # mosquito egg laying rate ("birth rate"); from paper
  phi_m = 1, # emergence rate; set equal to 1 to stabilize model dynamics
  p_m = 0.3, # probability of mosquito being infected by biting infectious human; set equal to p_h like paper
  # beta_m = 0.5 * 0.38, # human-to-mosquito transmission rate
  mu_m = 1/14 * 7, # adult mosquito mortality rate, set equal to egg laying rate
  mu_mL = 0.0 * 7, # larval mosquito death rate; initially set to 0
  p_EIP = 1/10 * 7 # progression rate out of mosquito E state; from paper
)

for (i in 1:length(model_parameters)) {
  julia_assign(names(model_parameters)[i], model_parameters[[i]])
}
#beta_h <- julia_eval("beta_h = c * p_h * ones(tmax);") # ((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/365))).+ 1);") #seasonal transmission- mosquito to human
beta_h <- julia_eval("beta_h = c * p_h * ((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/52))).+ 1);") #seasonal transmission- mosquito to human
beta_m <- julia_eval("beta_m = c * p_m;") # human-to-mosquito transmission rate

attach(model_parameters)
(m*median(beta_h)*beta_m*p_EIP)/(mu_m*(p_IP+bh)*(p_EIP+mu_m))
(m*c*p_h*c*p_m*p_EIP) / (mu_m*(p_IP+bh)*(p_EIP+mu_m))

m*c*p_m*1/p_IP*exp(-mu_m/p_EIP) * c*p_h*1/mu_m

c * p_m * 100 #lambdam
c * p_m * 100 * m*3255603 * 1/p_IP

# save all params into one vector to input into model function
julia_command("par = (; bh, beta_h, bm, phi_m, beta_m, mu_m, mu_mL, mu_h, p_IIP, p_IP, p_EIP, p_R, m)") 

######################################################
########### COMPARTMENTS
pop <- 3255603
julia_assign("pop", pop) # Initial total population of PR

compartment_names <- c('Sh_0', 'Eh_1', 'Ih_1', 'Rh_1', 'Sh_1', 
  'Eh_12', 'Eh_13', 'Eh_14', 'Ih_12', 'Ih_13', 'Ih_14', 
  'Rh_12', 'Rh_13', 'Rh_14', 
  'Eh_2', 'Ih_2', 'Rh_2', 'Sh_2', 
  'Eh_21', 'Eh_23', 'Eh_24', 'Ih_21', 'Ih_23', 'Ih_24', 'Rh_21', 'Rh_23', 'Rh_24', 
  'Eh_3', 'Ih_3', 'Rh_3', 'Sh_3', 
  'Eh_31', 'Eh_32', 'Eh_34', 'Ih_31', 'Ih_32', 'Ih_34', 'Rh_31', 'Rh_32', 'Rh_34', 
  'Eh_4', 'Ih_4', 'Rh_4', 'Sh_4', 
  'Eh_41', 'Eh_42', 'Eh_43', 'Ih_41', 'Ih_42', 'Ih_43', 
  'Rh_41', 'Rh_42', 'Rh_43', 
  'Ih_imp_1', 'Ih_imp_2', 'Ih_imp_3', 'Ih_imp_4',
  'Lm', 'Sm', 'Em1', 'Im1', 'Em2', 'Im2', 'Em3', 'Im3', 'Em4', 'Im4')
init_compartments <- matrix(0, nrow=length(compartment_names), ncol=nsims, 
  dimnames = list(compartment_names, paste0('V', 1:nsims)))

initial_s <- 0.5
initial_s2 <- rep(0.2/4, 4)
initial_r1 <- rep(0.2/4, 4)
initial_r2 <- rep(0.1/12, 12)
#init_h <- round(pop * t(rdirichlet(nsims, c(initial_s, initial_s2, initial_r1, initial_r2))))
# rdirichlet too variable, Sh_0 should always be high
init_h <- rmultinom(3, size=pop, prob=c(initial_s, initial_s2, initial_r1, initial_r2))
rownames(init_h) <- c('Sh_0', 'Sh_1', 'Sh_2', 'Sh_3', 'Sh_4', 
  "Rh_1", "Rh_2", "Rh_3", "Rh_4", "Rh_12", 'Rh_13', 'Rh_14', 'Rh_21', 'Rh_23', 'Rh_24', 
  'Rh_31', 'Rh_32', 'Rh_34', 'Rh_41', 'Rh_42', 'Rh_43')

init_compartments[rownames(init_h), ] <- init_h

# initialize mosquitos
init_compartments["Sm", ] <- pop * model_parameters$m

# use imported cases to initialize transmission
for (imp_compartment in c('Ih_imp_1', 'Ih_imp_2', 'Ih_imp_3', 'Ih_imp_4')) {
#  init_compartments[imp_compartment, ] <- rpois(nsims, 10)
  init_compartments[imp_compartment, ] <- rnbinom(nsims, mu=10, size = 0.1)
}

# transfer initial values to Julia
julia_assign("init_compartments", init_compartments)
julia_assign("input_compartment_names", compartment_names)
julia_command("input_compartment_names = Symbol.(input_compartment_names)")
julia_command("u0 = [];")
for (i in 1:nsims) {
  julia_command(paste0('push!(u0, NamedTuple{compartment_names}(init_compartments[:, ', i, ']))'))
}

######
# Vector of all compartments and sims
julia_command("outcomes0 = [];")
julia_command("for i=1:nsims
                push!(outcomes0, (; hpop = pop, mpop = pop * m, lpop = 0, 
                    hbirths = 0, hdeaths= 0, newcases_all_h = 0, newcases_all_m = 0, 
                    newcases_st1_h = 0, newcases_st2_h = 0, newcases_st3_h = 0, newcases_st4_h = 0,
                    p_infect_h1 = 0, p_infect_h2 = 0, p_infect_h3 = 0, p_infect_h4 = 0, 
                    p_infect_m1 = 0, p_infect_m2 = 0, p_infect_m3 = 0, p_infect_m4 = 0,
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

result <- julia_eval("run_4st_model_sims_randinit_imp!(nsims, tmax, u0, outcomes0, par)")

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
  plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Weeks", bty = "n")
  for (i in 1:nsims) {
    lines(1:tmax, df[, i], col = alpha(colors[1], 0.3))
  }
  lines(1:tmax, df$median, col = "black")
  
  x = 1:tmax
  
  for (i in 1:nsims) {
    y = df[, i]
    lw1 <- loess(y ~ x)
    # j <- order(df[, i])
    lines(x,lw1$fitted,col="red")
  }
}
#plot_all_humans()

# serotype 1 and 2 (humans)
plot_humans_by_seroptype <- function() {
  df1 <- as.data.frame(result$newcases_st1_h)
  df2 <- as.data.frame(result$newcases_st2_h)
  df3 <- as.data.frame(result$newcases_st3_h)
  df4 <- as.data.frame(result$newcases_st4_h)
  
  plot(NA, NA, xlim = c(0, tmax), ylim = c(10, max(df1, df2,df3,df4)), ylab = "New infections", xlab = "Weeks", bty = "n", log='y')
  for (i in 1:nsims) {
    lines(1:tmax, df1[, i], col = alpha(colors[1], 0.3))
  }
  # lines(1:tmax, df1$median, col = colors[2])
  for (i in 1:nsims) {
    lines(1:tmax, df2[, i], col = alpha(colors[2], 0.3))
  }
  # lines(1:tmax, df2$median, col = colors[3])
  for (i in 1:nsims) {
    lines(1:tmax, df3[, i], col = alpha(colors[3], 0.3))
  }
  for (i in 1:nsims) {
    lines(1:tmax, df4[, i], col = alpha(colors[4], 0.3))
  }
}
plot_humans_by_seroptype()


# df <- as.data.frame(I1dt)
# plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "output", xlab = "Days", bty = "n")
# for (i in 1:nsims) {
#   lines(1:tmax, df[, i], col = alpha(colors[2], 0.15))
# }


# serotype 2 (humans)
# df <- as.data.frame(newcases_st2_h)
# df$median <- apply(df, 1, median)
# plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Days", bty = "n")
# for (i in 1:nsims) {
#   lines(1:tmax, df[, i], col = alpha(colors[3], 0.15))
# }
# lines(1:tmax, df$median, col = colors[3])



############### FOI
# serotype 1 (humans)
plot_foi_serotype <- function(stype) {
  df <- as.data.frame(result[paste0('p_infect_h', stype)])
  df$median <- apply(df, 1, median)
  plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), 
    ylab = "p_inf_1", xlab = "Weeks", bty = "n")
  for (i in 1:nsims) {
    lines(1:tmax, df[, i], col = alpha(colors[4], 0.15))
  }
  lines(1:tmax, df$median, col = colors[4])
}
#plot_foi_serotype(1)

plot_populations <- function() {
  df <- as.data.frame(result$hpop)
  df$median <- apply(df, 1, median)
  plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "Human pop", xlab = "Weeks", bty = "n")
  for (i in 1:nsims) {
    lines(1:tmax, df[, i], col = alpha(colors[4], 0.15))
  }
  lines(1:tmax, df$median, col = colors[4])
  
  df <- as.data.frame(result$mpop)
  df$median <- apply(df, 1, median)
  plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "Mosquito pop", xlab = "Weeks", bty = "n")
  for (i in 1:nsims) {
    lines(1:tmax, df[, i], col = alpha(colors[4], 0.15))
  }
  lines(1:tmax, df$median, col = colors[4])
}
#plot_populations()

##### STACKED AREA PLOTS FOR COMPARTMENTS
source("data_stackedareaplot-MAJ.R")
data <- data_stackedareaplot2(list(
    Sh_0 = result$Sh_0, 
    Sh_secondary = result$Sh_1 + result$Sh_2 + result$Sh_3 + result$Sh_4,
    EI_primary = result$Eh_1 + result$Eh_2 + result$Eh_3 + result$Eh_4 + 
        result$Ih_1 + result$Ih_2 + result$Ih_3 + result$Ih_4,
    EI_secondary = result$Eh_12 + result$Eh_13 + result$Eh_14 + 
        result$Eh_21 + result$Eh_23 + result$Eh_24 + 
        result$Eh_31 + result$Eh_32 + result$Eh_34 +  
        result$Eh_41 + result$Eh_42 + result$Eh_43 +  
        result$Ih_12 + result$Ih_13 + result$Ih_14 + 
        result$Ih_21 + result$Ih_23 + result$Ih_24 + 
        result$Ih_31 + result$Ih_32 + result$Ih_34 +  
        result$Ih_41 + result$Ih_42 + result$Ih_43,
    R_primary = result$Rh_1 + result$Rh_2 + result$Rh_3 + result$Rh_4,
    R_secondary = result$Rh_12 + result$Rh_13 + result$Rh_14 + 
        result$Rh_21 + result$Rh_23 + result$Rh_24 + 
        result$Rh_31 + result$Rh_32 + result$Rh_34 +  
        result$Rh_41 + result$Rh_42 + result$Rh_43
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


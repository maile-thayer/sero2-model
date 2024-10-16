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
  p_IIP = min(1, 1/5 * 7), # progression rate out of human E state; from paper
  p_IP = min(1, 1/6 * 7), # progression rate out of human I state; from paper
  p_R = 1/365 * 7, # progression rate out of human cross-protective state; currently at 1 year; from paper
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

julia_command("par = (; bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R,m)") # save all params into one vector to input into model function

######################################################
########### COMPARTMENTS
pop <- 3255603
julia_assign("pop", pop) # Initial total population of PR
# initialize compartments
# susceptible and recovered estimated based on Sarah's catalytic model
# all E compartments empty to start
# start with 25 in each infectious human compartment

# init_comp_rand = Array{Float64}(undef, 9, nsims)
# 
# for (i in 1:nsims)(
#   
# )
# init_sims <- rdirichlet(nsims, c(0.2,0.025, 0.15, 0.05, 0.175,0.175, 0.025, 0.05, 0.15))

n_stype <- 4
initial_s <- 0.2#init_sims[,1]
initial_s2 <- (diff(c(0, sort(sample(seq(0.001, 0.999, 0.001), 3)), 1)))*(0.2) #rep(0.2/4, 4) #c(0.1, 0.1, 0.1, 0.1)#init_sims[,c(2:5)] #c(0.1, 0.1, 0.1, 0.1)#c(0.025, 0.15, 0.05, 0.175) #0.4
initial_r1 <- c(0.0, 0.0, 0.0, 0.0)
initial_r2 <- (diff(c(0, sort(sample(seq(0.001, 0.999, 0.001), 3)), 1)))*(0.6) #rep(0.6/4, 4)#c(0.1, 0.1, 0.1, 0.1)#init_sims[,c(6:9)]#c(0.1, 0.1, 0.1, 0.1)#c(0.175, 0.025, 0.05, 0.15) #0.4
(initial_s + sum(initial_s2) + sum(initial_r1) + sum(initial_r2))
init_sims <- rdirichlet(nsims, c(initial_s, initial_s2, initial_r2))

E_init <- matrix(0, nrow=(1+n_stype), ncol=n_stype)
I_init <- matrix(10, nrow=(1+n_stype), ncol=n_stype)
S_naive_init <- round(init_sims[,c(1)] * pop) - sum(I_init) - sum(E_init)
#Imp_init <- matrix(0, nrow=(1+n_stype), ncol=n_stype)
Imp_init <- matrix(0, nrow=1, ncol=n_stype)

# susceptible pop has to be > 0 --> keep re-drawing init_sims until this is true
if (sum(S_naive_init > 0) < nsims) {
  while (sum(S_naive_init > 0) < nsims) {
    init_sims <- rdirichlet(nsims, c(initial_s, initial_s2, initial_r2))
#    E_init <- matrix(0, nrow=(1+n_stype), ncol=n_stype)
#    I_init <- matrix(10, nrow=(1+n_stype), ncol=n_stype)
    S_naive_init <- round(init_sims[,c(1)] * pop) - sum(I_init) - sum(E_init)
  } 
} 

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

initial_s <- 0.2
initial_s2 <- rep(0.1/4, 4)
initial_r1 <- rep(0.3/4, 4)
initial_r2 <- rep(0.4/12, 12)
#init_h <- round(pop * t(rdirichlet(nsims, c(initial_s, initial_s2, initial_r1, initial_r2))))
# TODO: too variable, Sh_0 should always be high
init_h <- rmultinom(3, size=pop, prob=c(initial_s, initial_s2, initial_r1, initial_r2))
rownames(init_h) <- c('Sh_0', 'Sh_1', 'Sh_2', 'Sh_3', 'Sh_4', 
  "Rh_1", "Rh_2", "Rh_3", "Rh_4", "Rh_12", 'Rh_13', 'Rh_14', 'Rh_21', 'Rh_23', 'Rh_24', 
  'Rh_31', 'Rh_32', 'Rh_34', 'Rh_41', 'Rh_42', 'Rh_43')

init_compartments[rownames(init_h), ] <- init_h

# initialize mosquitos
init_compartments["Sm", ] <- pop * model_parameters$m

# use imported cases to initialize transmission
for (imp_compartment in c('Ih_imp_1', 'Ih_imp_2', 'Ih_imp_3', 'Ih_imp_4')) {
  init_compartments[imp_compartment, ] <- rpois(nsims, 10)
}

# transfer initial values to Julia
julia_assign("init_compartments", init_compartments)
julia_assign("input_compartment_names", compartment_names)
julia_command("input_compartment_names = Symbol.(input_compartment_names)")
julia_command("u0 = [];")
for (i in 1:nsims) {
  julia_command(paste0('push!(u0, NamedTuple{compartment_names}(init_compartments[:, ', i, ']))'))
}

# S_second_init <- round(init_sims[,c(2:5)] * pop)
# R_prim_init <- round(initial_r1 * pop)
# R_second_init <- array(NA,c(n_stype,n_stype,nsims))
# for (i in 1:nsims){
#   R_second_init[,,i] <- matrix(rep(round((init_sims[,c(6:9)][i,] * pop)/n_stype),n_stype), nrow=(n_stype), ncol=n_stype)
# }
# (sum(I_init) + sum(E_init) + S_naive_init + rowSums(S_second_init) + 
#     sum(R_prim_init) + apply(R_second_init,3,sum)) / pop
# 
# julia_assign('Sh_0', S_naive_init)
# julia_assign("Sh_1", S_second_init[,1])
# julia_assign("Sh_2", S_second_init[,2])
# julia_assign("Sh_3", S_second_init[,3])
# julia_assign("Sh_4", S_second_init[,4])
# julia_assign('Eh_1', E_init[1, 1])
# julia_assign('Eh_2', E_init[1, 2])
# julia_assign('Eh_3', E_init[1, 3])
# julia_assign('Eh_4', E_init[1, 4])
# julia_assign('Eh_12', E_init[2, 2])
# julia_assign('Eh_13', E_init[2, 3])
# julia_assign('Eh_14', E_init[2, 4])
# julia_assign('Eh_21', E_init[3, 1])
# julia_assign('Eh_23', E_init[3, 3])
# julia_assign('Eh_24', E_init[3, 4])
# julia_assign('Eh_31', E_init[4, 1])
# julia_assign('Eh_32', E_init[4, 2])
# julia_assign('Eh_34', E_init[4, 4])
# julia_assign('Eh_41', E_init[5, 1])
# julia_assign('Eh_42', E_init[5, 2])
# julia_assign('Eh_43', E_init[5, 3])
# julia_assign('Ih_1', I_init[1, 1])
# julia_assign('Ih_2', I_init[1, 2])
# julia_assign('Ih_3', I_init[1, 3])
# julia_assign('Ih_4', I_init[1, 4])
# julia_assign('Ih_12', I_init[2, 2])
# julia_assign('Ih_13', I_init[2, 3])
# julia_assign('Ih_14', I_init[2, 4])
# julia_assign('Ih_21', I_init[3, 1])
# julia_assign('Ih_23', I_init[3, 3])
# julia_assign('Ih_24', I_init[3, 4])
# julia_assign('Ih_31', I_init[4, 1])
# julia_assign('Ih_32', I_init[4, 2])
# julia_assign('Ih_34', I_init[4, 4])
# julia_assign('Ih_41', I_init[5, 1])
# julia_assign('Ih_42', I_init[5, 2])
# julia_assign('Ih_43', I_init[5, 3])
# julia_assign("Rh_1", R_prim_init[1])
# julia_assign("Rh_2", R_prim_init[2])
# julia_assign("Rh_3", R_prim_init[3])
# julia_assign("Rh_4", R_prim_init[4])
# julia_assign('Rh_12', R_second_init[1, 2, ])
# julia_assign('Rh_13', R_second_init[1, 3, ])
# julia_assign('Rh_14', R_second_init[1, 4, ])
# julia_assign('Rh_21', R_second_init[2, 1, ])
# julia_assign('Rh_23', R_second_init[2, 3, ])
# julia_assign('Rh_24', R_second_init[2, 4, ])
# julia_assign('Rh_31', R_second_init[3, 1, ])
# julia_assign('Rh_32', R_second_init[3, 2, ])
# julia_assign('Rh_34', R_second_init[3, 4, ])
# julia_assign('Rh_41', R_second_init[4, 1, ])
# julia_assign('Rh_42', R_second_init[4, 2, ])
# julia_assign('Rh_43', R_second_init[4, 3, ])
# 
# julia_assign('Ih_imp_1', Imp_init[1, 1])
# julia_assign('Ih_imp_2', Imp_init[1, 2])
# julia_assign('Ih_imp_3', Imp_init[1, 3])
# julia_assign('Ih_imp_4', Imp_init[1, 4])
# julia_assign('dIh_imp_12', Imp_init[2, 2])
# julia_assign('dIh_imp_13', Imp_init[2, 3])
# julia_assign('dIh_imp_14', Imp_init[2, 4])
# julia_assign('dIh_imp_21', Imp_init[3, 1])
# julia_assign('dIh_imp_23', Imp_init[3, 3])
# julia_assign('dIh_imp_24', Imp_init[3, 4])
# julia_assign('dIh_imp_31', Imp_init[4, 1])
# julia_assign('dIh_imp_32', Imp_init[4, 2])
# julia_assign('dIh_imp_34', Imp_init[4, 4])
# julia_assign('dIh_imp_41', Imp_init[5, 1])
# julia_assign('dIh_imp_42', Imp_init[5, 2])
# julia_assign('dIh_imp_43', Imp_init[5, 3])


# julia_command("dSh_0 = round(0.423*pop) ;")
# julia_command("dEh_1 = 0;")
# julia_command("dIh_1 = 25;")
# julia_command("dRh_1 = round((2/70)*0.227*pop);")
# julia_command("dEh_2 = 0;")
# julia_command("dIh_2 = 25;")
# julia_command("dRh_2 = round((2/70)*0.227*pop);")
# julia_command("dSh_1 = round((68/70)*0.227*pop);")
# julia_command("dEh_12 = 0;")
# julia_command("dIh_12 = 25;")
# julia_command("dRh_12 = round(0.061*pop);")
# julia_command("dSh_2 = round((68/70)*0.227*pop);")
# julia_command("dEh_21 = 0;")
# julia_command("dIh_21 = 25;")
# julia_command("dRh_21 = round(0.061*pop);")

# julia_command("Lm = 0;")
# # start with all mosquitoes (human pop * 2) in first S compartment
# julia_assign("Sm", pop*julia_eval('m'))
# julia_command("Em1 = 0;")
# julia_command("Im1 = 0;")
# julia_command("Em2 = 0;")
# julia_command("Im2 = 0;")
# julia_command("Em3 = 0;")
# julia_command("Im3 = 0;")
# julia_command("Em4 = 0;")
# julia_command("Im4 = 0;")
julia_command("sumdNm = sum(Sm+Em1+Im1+Em2+Im2+Em3+Im3+Em4+Im4) ;") # mosquito pop

######
# Vector of all compartments and sims
julia_command("outcomes0 = [];")
julia_command("for i=1:nsims
                push!(outcomes0, (; hpop = pop, mpop = sumdNm, lpop = 0, 
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

# test and run model 1 time (1 timestep)
# julia_command("x = dengue_2st!(x0,par,2)") #gives output
# x = julia_eval("dengue_2st!(x0,par,2)") #saves output as list

# set seed for reproducibility
# julia_command("Random.seed!(12345);")
# result <- julia_eval("run_4st_model_sims_randinit_imp!(nsims, tmax, x0, par)")
# 
# julia_source("/Users/michael/Library/CloudStorage/OneDrive-Personal/Documents/Projects/learn julia/test_named_output.jl")
# result <- julia_eval("test_named_output!(nsims, tmax, x0, par)")

result <- julia_eval("run_4st_model_sims_randinit_imp!(nsims, tmax, u0, outcomes0, par)")
head(result$newcases_all_h)

###############################################################
###############################################################
##### 4. Plot output
###############################################################
###############################################################
colors <- brewer.pal(8, name = "Dark2")

############### NEWCASES
# all (humans)
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


# serotype 1 and 2 (humans)
df1 <- as.data.frame(result$newcases_st1_h)
df2 <- as.data.frame(result$newcases_st2_h)
df3 <- as.data.frame(result$newcases_st3_h)
df4 <- as.data.frame(result$newcases_st4_h)
# df1$median <- apply(df1, 1, median)
# df2$median <- apply(df2, 1, median)
# df1 <- as.data.frame(I1dt + I21dt + I31dt + I41dt)
# df2 <- as.data.frame(I2dt + I12dt + I32dt + I42dt)
# df3 <- as.data.frame(I3dt + I13dt + I23dt + I43dt)
# df4 <- as.data.frame(I4dt + I14dt + I24dt + I34dt)
# df1 <- as.data.frame(I1mdt)
# df2 <- as.data.frame(I2mdt)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df1, df2,df3,df4)), ylab = "New infections", xlab = "Weeks", bty = "n")
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
df <- as.data.frame(result$p_infect_h1)
df$median <- apply(df, 1, median)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), 
  ylab = "p_inf_1", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df[, i], col = alpha(colors[4], 0.15))
}
lines(1:tmax, df$median, col = colors[4])

# serotype 2 (humans)
df <- as.data.frame(result$p_infect_h2)
df$median <- apply(df, 1, median)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "p_inf_2", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df[, i], col = alpha(colors[5], 0.15))
}
lines(1:tmax, df$median, col = colors[5])


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


##### STACKED AREA PLOTS FOR COMPARTMENTS

S.col <- brewer.pal(9, "GnBu")[7:9]
E.col <- brewer.pal(9, "Purples")[6:9]
I.col <- brewer.pal(9, "YlOrRd")[6:9]
R.col <- brewer.pal(9, "BuGn")[6:9]

source("data_stackedareaplot.R")


data <- data_stackedareaplot(result$Sh_0, result$Sh_1, result$Sh_2,
  result$Eh_1, result$Eh_2, result$Eh_12, result$Eh_21,
  result$Ih_1, result$Ih_2, result$Ih_12, result$Ih_21,
  result$Rh_1, result$Rh_2, result$Rh_12, result$Rh_21,
  nsims = nsims, tmax = tmax
)

# Population distribution of each compartment over timeframe
ggplot(data, aes(x = time, y = value_perc, fill = group)) +
  geom_area() +
  theme_classic() +
  ggtitle("4-serotype model") +
  scale_fill_manual(values = c(R.col, S.col, E.col, I.col))


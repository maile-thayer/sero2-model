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
julia <- julia_setup(JULIA_HOME = "C:/Users/ruu6/AppData/Local/Programs/Julia-1.7.2/bin")

############################################
# setwd("~/GitHub/sero2-model/") #set working directory to github project repo if you haven't already
# julia_install_package("Random")
# julia_install_package_if_needed("Random")
julia_library("Random, Distributions, DataFrames, LinearAlgebra, CSV, JLD2, Dates")

julia_source("1_dengue_4st_imp_diffp.jl") # load main model function
julia_source("2_run_4st_model_sims_randinit_imp_diffp.jl") # load function to run main model

julia_exists("dengue_4st_imp_diffp!") # test whether functions loaded and exist :)
julia_exists("multinom_samp")
julia_exists("run_4st_model_sims_randinit_imp_diffp!")

##############################################################
###############################################################
##### 2. Initialize parameters and populations/compartments
###############################################################
###############################################################

##### PARAMETERS
# set parameters using Wearing et al paper https://www.pnas.org/doi/full/10.1073/pnas.0602960103#:~:text=Rapid%20and%20unplanned%20urbanization%2C%20increased,contributed%20to%20dramatic%20dengue%20resurgence.


# set timeframe
tmax <- julia_eval("tmax = 52*1000;")
julia_command("tspan = collect(0.0:(tmax+52));")

# set # simulations to run-- doing it in R
# nsims=2
nsims <- julia_eval("nsims = 20;")


detach(model_parameters)
######## HUMAN PARAMETERS ########
model_parameters <- list(
  # human
  bh = 1/(60*365) * 7, # birthrate; from paper
  m = 2, # initial ratio of mosquitoes to humans; from literature
  c = 0.5 * 7, # contact rate (biting rate); from paper
  p_h1 = 0.1, # probability of a human being infected by an infected mosquito bite; paper 0.38
  p_h2 = 0.09,
  p_h3 = 0.11,
  p_h4 = 0.1,
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

# p_h <- julia_assign("p_h", "[p_h1 p_h2 p_h3 p_h4]")
#beta_h <- julia_eval("beta_h = c * p_h * ones(tmax);") # ((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/365))).+ 1);") #seasonal transmission- mosquito to human
# beta_h <- julia_eval("beta_h = c * p_h * ((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/52))).+ 1);") #seasonal transmission- mosquito to human
beta_h <- julia_eval("beta_h = c * [p_h1 p_h2 p_h3 p_h4].* ((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/52))).+ 1);")
beta_m <- julia_eval("beta_m = c * p_m;") # human-to-mosquito transmission rate

attach(model_parameters)
(m*median(beta_h)*beta_m*p_EIP)/(mu_m*(p_IP+bh)*(p_EIP+mu_m))
# (m*c*p_h*c*p_m*p_EIP) / (mu_m*(p_IP+bh)*(p_EIP+mu_m))

(m*c*c(p_h1,p_h2,p_h3,p_h4)*c*p_m*p_EIP) / (mu_m*(p_IP+bh)*(p_EIP+mu_m))

m*c*p_m*1/p_IP*exp(-mu_m/p_EIP) * c*c(p_h1,p_h2,p_h3,p_h4)*1/mu_m

c * p_m * 100 #lambdam
c * p_m * 100 * m*3255603 * 1/p_IP

julia_command("par = (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R,m);") # save all params into one vector to input into model function

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
Imp_init <- matrix(0, nrow=(1+n_stype), ncol=n_stype)

# susceptible pop has to be > 0 --> keep re-drawing init_sims until this is true
if (sum(S_naive_init > 0) < nsims) {
  while (sum(S_naive_init > 0) < nsims) {
    init_sims <- rdirichlet(nsims, c(initial_s, initial_s2, initial_r2))
    E_init <- matrix(0, nrow=(1+n_stype), ncol=n_stype)
    I_init <- matrix(10, nrow=(1+n_stype), ncol=n_stype)
    S_naive_init <- round(init_sims[,c(1)] * pop) - sum(I_init) - sum(E_init)
  } 
} 

S_second_init <- round(init_sims[,c(2:5)] * pop)
R_prim_init <- round(initial_r1 * pop)
R_second_init <- array(NA,c(n_stype,n_stype,nsims))
for (i in 1:nsims){
  R_second_init[,,i] <- matrix(rep(round((init_sims[,c(6:9)][i,] * pop)/n_stype),n_stype), nrow=(n_stype), ncol=n_stype)
}
(sum(I_init) + sum(E_init) + S_naive_init + rowSums(S_second_init) + 
    sum(R_prim_init) + apply(R_second_init,3,sum)) / pop

julia_assign('u0Sh_0', S_naive_init)
julia_assign("u0Sh_1", S_second_init[,1])
julia_assign("u0Sh_2", S_second_init[,2])
julia_assign("u0Sh_3", S_second_init[,3])
julia_assign("u0Sh_4", S_second_init[,4])
julia_assign('u0Eh_1', E_init[1, 1])
julia_assign('u0Eh_2', E_init[1, 2])
julia_assign('u0Eh_3', E_init[1, 3])
julia_assign('u0Eh_4', E_init[1, 4])
julia_assign('u0Eh_12', E_init[2, 2])
julia_assign('u0Eh_13', E_init[2, 3])
julia_assign('u0Eh_14', E_init[2, 4])
julia_assign('u0Eh_21', E_init[3, 1])
julia_assign('u0Eh_23', E_init[3, 3])
julia_assign('u0Eh_24', E_init[3, 4])
julia_assign('u0Eh_31', E_init[4, 1])
julia_assign('u0Eh_32', E_init[4, 2])
julia_assign('u0Eh_34', E_init[4, 4])
julia_assign('u0Eh_41', E_init[5, 1])
julia_assign('u0Eh_42', E_init[5, 2])
julia_assign('u0Eh_43', E_init[5, 3])
julia_assign('u0Ih_1', I_init[1, 1])
julia_assign('u0Ih_2', I_init[1, 2])
julia_assign('u0Ih_3', I_init[1, 3])
julia_assign('u0Ih_4', I_init[1, 4])
julia_assign('u0Ih_12', I_init[2, 2])
julia_assign('u0Ih_13', I_init[2, 3])
julia_assign('u0Ih_14', I_init[2, 4])
julia_assign('u0Ih_21', I_init[3, 1])
julia_assign('u0Ih_23', I_init[3, 3])
julia_assign('u0Ih_24', I_init[3, 4])
julia_assign('u0Ih_31', I_init[4, 1])
julia_assign('u0Ih_32', I_init[4, 2])
julia_assign('u0Ih_34', I_init[4, 4])
julia_assign('u0Ih_41', I_init[5, 1])
julia_assign('u0Ih_42', I_init[5, 2])
julia_assign('u0Ih_43', I_init[5, 3])
julia_assign("u0Rh_1", R_prim_init[1])
julia_assign("u0Rh_2", R_prim_init[2])
julia_assign("u0Rh_3", R_prim_init[3])
julia_assign("u0Rh_4", R_prim_init[4])
julia_assign('u0Rh_12', R_second_init[1, 2, ])
julia_assign('u0Rh_13', R_second_init[1, 3, ])
julia_assign('u0Rh_14', R_second_init[1, 4, ])
julia_assign('u0Rh_21', R_second_init[2, 1, ])
julia_assign('u0Rh_23', R_second_init[2, 3, ])
julia_assign('u0Rh_24', R_second_init[2, 4, ])
julia_assign('u0Rh_31', R_second_init[3, 1, ])
julia_assign('u0Rh_32', R_second_init[3, 2, ])
julia_assign('u0Rh_34', R_second_init[3, 4, ])
julia_assign('u0Rh_41', R_second_init[4, 1, ])
julia_assign('u0Rh_42', R_second_init[4, 2, ])
julia_assign('u0Rh_43', R_second_init[4, 3, ])

julia_assign('u0Ih_imp_1', Imp_init[1, 1])
julia_assign('u0Ih_imp_2', Imp_init[1, 2])
julia_assign('u0Ih_imp_3', Imp_init[1, 3])
julia_assign('u0Ih_imp_4', Imp_init[1, 4])
julia_assign('u0Ih_imp_12', Imp_init[2, 2])
julia_assign('u0Ih_imp_13', Imp_init[2, 3])
julia_assign('u0Ih_imp_14', Imp_init[2, 4])
julia_assign('u0Ih_imp_21', Imp_init[3, 1])
julia_assign('u0Ih_imp_23', Imp_init[3, 3])
julia_assign('u0Ih_imp_24', Imp_init[3, 4])
julia_assign('u0Ih_imp_31', Imp_init[4, 1])
julia_assign('u0Ih_imp_32', Imp_init[4, 2])
julia_assign('u0Ih_imp_34', Imp_init[4, 4])
julia_assign('u0Ih_imp_41', Imp_init[5, 1])
julia_assign('u0Ih_imp_42', Imp_init[5, 2])
julia_assign('u0Ih_imp_43', Imp_init[5, 3])


# julia_command("u0Sh_0 = round(0.423*pop) ;")
# julia_command("u0Eh_1 = 0;")
# julia_command("u0Ih_1 = 25;")
# julia_command("u0Rh_1 = round((2/70)*0.227*pop);")
# julia_command("u0Eh_2 = 0;")
# julia_command("u0Ih_2 = 25;")
# julia_command("u0Rh_2 = round((2/70)*0.227*pop);")
# julia_command("u0Sh_1 = round((68/70)*0.227*pop);")
# julia_command("u0Eh_12 = 0;")
# julia_command("u0Ih_12 = 25;")
# julia_command("u0Rh_12 = round(0.061*pop);")
# julia_command("u0Sh_2 = round((68/70)*0.227*pop);")
# julia_command("u0Eh_21 = 0;")
# julia_command("u0Ih_21 = 25;")
# julia_command("u0Rh_21 = round(0.061*pop);")

julia_command("u0Lm = 0;")
# start with all mosquitoes (human pop * 2) in first S compartment
julia_assign("u0Sm", pop*julia_eval('m'))
julia_command("u0Em1 = 0;")
julia_command("u0Im1 = 0;")
julia_command("u0Em2 = 0;")
julia_command("u0Im2 = 0;")
julia_command("u0Em3 = 0;")
julia_command("u0Im3 = 0;")
julia_command("u0Em4 = 0;")
julia_command("u0Im4 = 0;")
julia_command("sumu0Nm = sum(u0Sm+u0Em1+u0Im1+u0Em2+u0Im2+u0Em3+u0Im3+u0Em4+u0Im4) ;") # mosquito pop

######
# Vector of all compartments and sims
julia_command("u0all = [];")
julia_command("x0 = [];")
julia_command("for j=1:nsims
                push!(u0all,[u0Sh_0[j],u0Eh_1,u0Ih_1,u0Rh_1,u0Sh_1[j],u0Eh_12,u0Eh_13,u0Eh_14,u0Ih_12,u0Ih_13,u0Ih_14,u0Rh_12[j],u0Rh_13[j],u0Rh_14[j], 
                                 u0Eh_2,u0Ih_2,u0Rh_2,u0Sh_2[j],u0Eh_21,u0Eh_23,u0Eh_24,u0Ih_21,u0Ih_23,u0Ih_24,u0Rh_21[j],u0Rh_23[j],u0Rh_24[j], 
                                 u0Eh_3,u0Ih_3,u0Rh_3,u0Sh_3[j],u0Eh_31,u0Eh_32,u0Eh_34,u0Ih_31,u0Ih_32,u0Ih_34,u0Rh_31[j],u0Rh_32[j],u0Rh_34[j], 
                                 u0Eh_4,u0Ih_4,u0Rh_4,u0Sh_4[j],u0Eh_41,u0Eh_42,u0Eh_43,u0Ih_41,u0Ih_42,u0Ih_43,u0Rh_41[j],u0Rh_42[j],u0Rh_43[j], 
                                 u0Ih_imp_1,u0Ih_imp_2,u0Ih_imp_3,u0Ih_imp_4,
                                 u0Ih_imp_12,u0Ih_imp_13,u0Ih_imp_14,
                                 u0Ih_imp_21,u0Ih_imp_23,u0Ih_imp_24,
                                 u0Ih_imp_31,u0Ih_imp_32,u0Ih_imp_34,
                                 u0Ih_imp_41,u0Ih_imp_42,u0Ih_imp_43, 
                                 u0Lm,u0Sm,u0Em1,u0Im1,u0Em2,u0Im2,u0Em3,u0Im3,u0Em4,u0Im4])
                push!(x0,[u0all[j],0,0,pop,sumu0Nm,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
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
julia_command("Random.seed!(12345);")
result <- julia_eval("run_4st_model_sims_randinit_imp_diffp!(nsims, tmax, x0, par)")
# saved result is a 'JuliaObject' list of 34 nsims*tmax matrices

# CASES
newcases_all_h <- result[[1]] # all new human infections/cases from time 1:tmax
newcases_all_m <- result[[2]] # all new mosquito infections/cases from time 1:tmax
newcases_st1_h <- result[[3]] # new serotype-1 human infections/cases from time 1:tmax
newcases_st2_h <- result[[4]] # new serotype-2 human infections/cases from time 1:tmax
newcases_st3_h <- result[[5]] # new serotype-3 human infections/cases from time 1:tmax
newcases_st4_h <- result[[6]] # new serotype-4 human infections/cases from time 1:tmax
# POPULATIONS, BIRTHS, DEATHS
hpop <- result[[7]] # human population from time 1:tmax
mpop <- result[[8]] # adult mosquito population from time 1:tmax
lpop <- result[[9]] # larval mosquito population from time 1:tmax
hbirths <- result[[10]] # human births from time 1:tmax
hdeaths <- result[[11]] # human deaths from time 1:tmax

# SEROTYPE-SPECIFIC FORCES OF INFECTION
lambda_h1 <- result[[12]] # mosquito-to-human serotype-1 FOI
lambda_h2 <- result[[13]] # mosquito-to-human serotype-2 FOI
lambda_h2 <- result[[14]] # mosquito-to-human serotype-3 FOI
lambda_h2 <- result[[15]] # mosquito-to-human serotype-4 FOI
lambda_m1 <- result[[16]] # human-to-mosquito serotype-1 FOI
lambda_m2 <- result[[17]] # human-to-mosquito serotype-2 FOI
lambda_m3 <- result[[18]] # human-to-mosquito serotype-3 FOI
lambda_m4 <- result[[19]] # human-to-mosquito serotype-4 FOI

# HUMAN COMPARTMENTS
S0dt <- result[[20]] # Susceptible to any serotype (0 prior infections)
E1dt <- result[[21]] # Incubation compartment (IIP) after being infected by st1 for first time (only st1; no others)
E2dt <- result[[22]] # Incubation compartment (IIP) after being infected by st2 for first time (only st2; no others)
E3dt <- result[[23]] # Incubation compartment (IIP) after being infected by st3 for first time (only st3; no others)
E4dt <- result[[24]] # Incubation compartment (IIP) after being infected by st4 for first time (only st4; no others)
I1dt <- result[[25]] # Infectious compartment (IP) for st1 (only st1; no others)
I2dt <- result[[26]] # Infectious compartment (IP) for st2 (only st2; no others)
I3dt <- result[[27]] # Infectious compartment (IP) for st3 (only st3; no others)
I4dt <- result[[28]] # Infectious compartment (IP) for st4 (only st4; no others)
R1dt <- result[[29]] # Recovered and temporarily immune (cross-protection) compartment for st1
R2dt <- result[[30]] # Recovered and temporarily immune (cross-protection) compartment for st2
R3dt <- result[[31]] # Recovered and temporarily immune (cross-protection) compartment for st3
R4dt <- result[[32]] # Recovered and temporarily immune (cross-protection) compartment for st4
S1dt <- result[[33]] # Susceptible to st2,3,4 (but has been previously infected with st1)
S2dt <- result[[34]] # Susceptible to st1,3,4 (but has been previously infected with st2)
S3dt <- result[[35]] # Susceptible to st1,2,4 (but has been previously infected with st3)
S4dt <- result[[36]] # Susceptible to st1,2,3 (but has been previously infected with st4)
E12dt <- result[[37]] # Incubation compartment (IIP) after being infected by st2 (previously infected by st1)
E13dt <- result[[38]] # Incubation compartment (IIP) after being infected by st3 (previously infected by st1)
E14dt <- result[[39]] # Incubation compartment (IIP) after being infected by st4 (previously infected by st1)
I12dt <- result[[40]] # Infectious compartment (IP) for st2 (previously infected by st1)
I13dt <- result[[41]] # Infectious compartment (IP) for st3 (previously infected by st1)
I14dt <- result[[42]] # Infectious compartment (IP) for st4 (previously infected by st1)
R12dt <- result[[43]] # Recovered and immune from infection with serotypes 1,2
R13dt <- result[[44]] # Recovered and immune from infection with serotypes 1,3
R14dt <- result[[45]] # Recovered and immune from infection with serotypes 1,4
E21dt <- result[[46]] # Incubation compartment (IIP) after being infected by st1 (previously infected by st2)
E23dt <- result[[47]] # Incubation compartment (IIP) after being infected by st3 (previously infected by st2)
E24dt <- result[[48]] # Incubation compartment (IIP) after being infected by st4 (previously infected by st2)
I21dt <- result[[49]] # Infectious compartment (IP) for st1 (previously infected by st2)
I23dt <- result[[50]] # Infectious compartment (IP) for st3 (previously infected by st2)
I24dt <- result[[51]] # Infectious compartment (IP) for st4 (previously infected by st2)
R21dt <- result[[52]] # Recovered and immune from infection with serotypes 2,1
R23dt <- result[[53]] # Recovered and immune from infection with serotypes 2,3
R24dt <- result[[54]] # Recovered and immune from infection with serotypes 2,4
E31dt <- result[[55]] # Incubation compartment (IIP) after being infected by st1 (previously infected by st3)
E32dt <- result[[56]] # Incubation compartment (IIP) after being infected by st2 (previously infected by st3)
E34dt <- result[[57]] # Incubation compartment (IIP) after being infected by st4 (previously infected by st3)
I31dt <- result[[58]] # Infectious compartment (IP) for st1 (previously infected by st3)
I32dt <- result[[59]] # Infectious compartment (IP) for st2 (previously infected by st3)
I34dt <- result[[60]] # Infectious compartment (IP) for st4 (previously infected by st3)
R31dt <- result[[61]] # Recovered and immune from infection with serotypes 3,1
R32dt <- result[[62]] # Recovered and immune from infection with serotypes 3,2
R34dt <- result[[63]] # Recovered and immune from infection with serotypes 3,4
E41dt <- result[[64]] # Incubation compartment (IIP) after being infected by st1 (previously infected by st4)
E42dt <- result[[65]] # Incubation compartment (IIP) after being infected by st2 (previously infected by st4)
E43dt <- result[[66]] # Incubation compartment (IIP) after being infected by st3 (previously infected by st4)
I41dt <- result[[67]] # Infectious compartment (IP) for st1 (previously infected by st4)
I42dt <- result[[68]] # Infectious compartment (IP) for st1 (previously infected by st4)
I43dt <- result[[69]] # Infectious compartment (IP) for st1 (previously infected by st4)
R41dt <- result[[70]] # Recovered and immune from infection with serotypes 4,1
R42dt <- result[[71]] # Recovered and immune from infection with serotypes 4,2
R43dt <- result[[72]] # Recovered and immune from infection with serotypes 4,3

# MOSQUITO COMPARTMENTS -- can only be infected once
Lmdt <- result[[73]] # Development (Larval) compartment (no infections)
Smdt <- result[[74]] # Susceptible to any serotype
E1mdt <- result[[75]] # Incubation compartment (EIP) after being infected by st1
E2mdt <- result[[76]] # Incubation compartment (EIP) after being infected by st2
E3mdt <- result[[77]] # Incubation compartment (EIP) after being infected by st3
E4mdt <- result[[78]] # Incubation compartment (EIP) after being infected by st4
I1mdt <- result[[79]] # Infectious compartment after being infected by st1
I2mdt <- result[[80]] # Infectious compartment after being infected by st2
I3mdt <- result[[81]] # Infectious compartment after being infected by st3
I4mdt <- result[[82]] # Infectious compartment after being infected by st4
Impdt <- result[[83]] # 
am_pop <- Smdt + E1mdt + E2mdt + E3mdt + E4mdt + I1mdt + I2mdt + I3mdt + I4mdt

###############################################################
###############################################################
##### 4. Plot output
###############################################################
###############################################################
colors <- brewer.pal(8, name = "Dark2")


############### NEWCASES
# all (humans)
df <- as.data.frame(newcases_all_h)
df$median <- apply(df, 1, median)

# i = 1
# plot(NA, NA, xlim = c(4000, 5000), ylim = c(0, 20), ylab = "New infections", xlab = "Weeks", bty = "n")
# lines(1:tmax, df[, 1], col = colors[1])
# lines(1:tmax, df[, 2], col = colors[2])
# lines(1:tmax, df[, 3], col = colors[3])
j = j+1
plot(NA, NA, xlim = c(5000+(j*1000), 5000+(j*1000)+1000), ylim = c(0, 30), ylab = "New infections", xlab = "Weeks", bty = "n")
lines(1:tmax, df[, 1], col = colors[1])



plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df[, i], col = alpha(colors[1], 0.3))
}
lines(1:tmax, df$median, col = "black", lwd=2)

x = 1:tmax

for (i in 1:nsims) {
  y = df[, i]
  lw1 <- loess(y ~ x)
  # j <- order(df[, i])
  lines(x,lw1$fitted,col="red")
}


# serotype 1 and 2 (humans)
df1 <- as.data.frame(newcases_st1_h)
df2 <- as.data.frame(newcases_st2_h)
df3 <- as.data.frame(newcases_st3_h)
df4 <- as.data.frame(newcases_st4_h)
df1$median <- apply(df1, 1, median)
df2$median <- apply(df2, 1, median)
df3$median <- apply(df3, 1, median)
df4$median <- apply(df4, 1, median)
# df1 <- as.data.frame(I1dt + I21dt + I31dt + I41dt)
# df2 <- as.data.frame(I2dt + I12dt + I32dt + I42dt)
# df3 <- as.data.frame(I3dt + I13dt + I23dt + I43dt)
# df4 <- as.data.frame(I4dt + I14dt + I24dt + I34dt)
# df1 <- as.data.frame(I1mdt)
# df2 <- as.data.frame(I2mdt)
# plot(NA, NA, xlim = c(tmax-1000, tmax), ylim = c(0, 20), ylab = "New infections", xlab = "Weeks", bty = "n")
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, 5), ylab = "New infections", xlab = "Weeks", bty = "n")

plot(NA, NA, xlim = c(1000, 2000), ylim = c(0, 5), ylab = "New infections", xlab = "Weeks", bty = "n")

plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df1, df2,df3,df4)), ylab = "New infections", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df1[, i], col = alpha(colors[1], 0.1))
}
# lines(1:tmax, df1$median, col = colors[2])
for (i in 1:nsims) {
  lines(1:tmax, df2[, i], col = alpha(colors[2], 0.1))
}
# lines(1:tmax, df2$median, col = colors[3])
for (i in 1:nsims) {
  lines(1:tmax, df3[, i], col = alpha(colors[3], 0.1))
}
for (i in 1:nsims) {
  lines(1:tmax, df4[, i], col = alpha(colors[4], 0.1))
}
lines(1:tmax, df1$median, col = colors[1])
lines(1:tmax, df2$median, col = colors[2])
lines(1:tmax, df3$median, col = colors[3])
lines(1:tmax, df4$median, col = colors[4])

lines(1:tmax, df1$median, col = "black")
lines(1:tmax, df2$median, col = "black")
lines(1:tmax, df3$median, col = "black")
lines(1:tmax, df4$median, col = "black")

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
df <- as.data.frame(lambda_h1)
df$median <- apply(df, 1, median)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df[, i], col = alpha(colors[4], 0.15))
}
lines(1:tmax, df$median, col = colors[4])

# serotype 2 (humans)
df <- as.data.frame(lambda_h2)
df$median <- apply(df, 1, median)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df[, i], col = alpha(colors[5], 0.15))
}
lines(1:tmax, df$median, col = colors[5])


X <- hpop
df <- as.data.frame(X)
df$median <- apply(df, 1, median)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Weeks", bty = "n")
for (i in 1:nsims) {
  lines(1:tmax, df[, i], col = alpha(colors[4], 0.15))
}
lines(1:tmax, df$median, col = colors[4])

X <- mpop
df <- as.data.frame(X)
df$median <- apply(df, 1, median)
plot(NA, NA, xlim = c(0, tmax), ylim = c(0, max(df)), ylab = "New infections", xlab = "Weeks", bty = "n")
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

data <- data_stackedareaplot(S0dt, S1dt, S2dt,
  E1dt, E2dt, E12dt, E21dt,
  I1dt, I2dt, I12dt, I21dt,
  R1dt, R2dt, R12dt, R21dt,
  nsims = nsims, tmax = tmax
)

# Population distribution of each compartment over timeframe
ggplot(data, aes(x = time, y = value_perc, fill = group)) +
  geom_area() +
  theme_classic() +
  ggtitle("2-serotype model") +
  scale_fill_manual(values = c(R.col, S.col, E.col, I.col))


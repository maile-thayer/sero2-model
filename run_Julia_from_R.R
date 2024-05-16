
###############################################################
###############################################################
##### 1. Install JuliaCall in R, load packages and environment
###############################################################
###############################################################

# install.packages("JuliaCall") 

library(JuliaCall)

# install_julia() #DONT RUN UNLESS YOU DON'T HAVE JULIA INSTALLED
julia <- julia_setup(JULIA_HOME = "C:/Users/ruu6/AppData/Local/Programs/Julia-1.7.2/bin")

############################################
# need "julia_..." in front of command to invoke program
# put commands in parentheses ""
#test out sqrt(2)
julia_eval("sqrt(2)")
#assign
julia_command("a = sqrt(2);")
julia_eval("a")

# julia documentation: https://docs.julialang.org/en/v1/ 


############################################
# setwd("~/GitHub/sero2-model/") #set to github project repo dir if you haven't already
# julia_install_package("Random")
# julia_install_package_if_needed("Random")
julia_library("Random")
julia_library("Distributions")
julia_library("DataFrames")
julia_library("LinearAlgebra")
julia_library("CSV")
julia_library("JLD2")
julia_library("Dates")

julia_source("1_dengue_2st.jl") #load main model function
julia_source("2_run_model_sims.jl") #load function to raun main model

julia_exists("dengue_2st!")
julia_exists("dengue_2st!")

# using Random,Distributions,DataFrames,LinearAlgebra,CSV,JLD2,Dates


###############################################################
###############################################################
##### 2. Initialize parameters and initial populations
###############################################################
###############################################################
# set parameters using Wearing et al paper https://www.pnas.org/doi/full/10.1073/pnas.0602960103#:~:text=Rapid%20and%20unplanned%20urbanization%2C%20increased,contributed%20to%20dramatic%20dengue%20resurgence. 

# set seed for reproducibility
julia_command("Random.seed!(123);")

# set timeframe
julia_command("tmax = 3650;")
julia_command("tspan = collect(0.0:(tmax+500));")

#set # simulations to run
julia_command("nsims = 10;")

######## HUMAN PARAMETERS ########
julia_command("bh = 1/(60*365);") #birthrate; from paper
julia_command("m = 2;") #initial ratio of mosquitoes to humans; from literature
julia_command("c = 0.5;") #contact rate (biting rate); from paper
julia_command("p_h = 0.44;") #probability of a human being infected by an infected mosquito bite; paper 0.38
julia_command("beta_h = c*p_h*((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/365))).+ 1);") #seasonal transmission- mosquito to human
julia_command("mu_h = bh;") #human death rate; set equal to birthrate
julia_command("p_IIP = 1/5;") #progression rate out of human E state; from paper

julia_command("p_IP = 1/6;")#progression rate out of human I state; from paper

julia_command("p_R = 1/365;") #progression rate out of human cross-protective state; currently at 1 year; from paper


######## MOSQUITO PARAMS ########
julia_command("bm = 1/14;") #mosquito egg laying rate ("birth rate"); from paper
julia_command("phi_m = 1;") #emergence rate; set equal to 1 to stabilize model dynamics
julia_command("p_m = p_h;") #probability of mosquito being infected by biting infectious human; set equal to p_h like paper
julia_command("beta_m = c*p_m;") # human-to-mosquito transmission rate
julia_command("mu_m = bm;") #adult mosquito mortality rate, set equal to egg laying rate
julia_command("mu_mL = 0.0;") #larval mosquito death rate; initially set to 0
julia_command("p_EIP = 1/10;") #progression rate out of mosquito E state; from paper



julia_command("R0 = (m*median(beta_h)*beta_m*p_EIP)/(mu_m*(p_IP+bh)*(p_EIP+mu_m));") #check R0
julia_eval("R0")


julia_command("par = (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R);") #save all params into one vector to input into model function

######################################################
# Initialize compartment pop
julia_command("pop = 3255603;")#Initial total population of PR
#initialize compartments
#susceptible and recovered estimated based on Sarah's catalytic model
#all E compartments empty to start
#start with 25 in each infectious human compartment
julia_command("u0Sh_0 = round(0.423*pop) ;")
julia_command("u0Eh_1 = 0;")
julia_command("u0Ih_1 = 25;")
julia_command("u0Rh_1 = round((2/70)*0.227*pop);")
julia_command("u0Eh_2 = 0;")
julia_command("u0Ih_2 = 25;")
julia_command("u0Rh_2 = round((2/70)*0.227*pop);")
julia_command("u0Sh_1 = round((68/70)*0.227*pop);")
julia_command("u0Eh_12 = 0;")
julia_command("u0Ih_12 = 25;")
julia_command("u0Rh_12 = round(0.061*pop);")
julia_command("u0Sh_2 = round((68/70)*0.227*pop);")
julia_command("u0Eh_21 = 0;")
julia_command("u0Ih_21 = 25;")
julia_command("u0Rh_21 = round(0.061*pop);")
julia_command("u0Lm = 0;")
#start with all mosquitoes (human pop * 2) in first S compartment
julia_command("u0Sm = pop*m;")
julia_command("u0Em1 = 0;")
julia_command("u0Im1 = 0;")
julia_command("u0Em2 = 0;")
julia_command("u0Im2 = 0;")

julia_command("sumu0Nm = sum(u0Sm+u0Em1+u0Im1+u0Em2+u0Im2) ;") #mosquito pop

######
# Vector of all compartments
julia_command("u0all = ([u0Sh_0, u0Eh_1, u0Eh_2, u0Ih_1, u0Ih_2, u0Rh_1, u0Rh_2, u0Sh_1, u0Sh_2, u0Eh_12, u0Eh_21, u0Ih_12, u0Ih_21, u0Rh_12, u0Rh_21, u0Lm, u0Sm, u0Em1, u0Em2, u0Im1, u0Im2]);")
# initial value for all later stored compartments (see below)
julia_command("x0 = [u0all,0,0,pop,sumu0Nm,0,0,0,0,0,0,0,0,0,0,0,0,0,0];")
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
# 10. mosquito force of infection st1
# 11. mosquito force of infection st2
# 12. all new FIRST infections
# X 13. aging (by age) ---- not using in this model
# X 14. deaths by age ---- not using in this model
# X 15. importations ---- not using in this model
# 16. new secondary human cases
# 17. total larval population
# 18. new human serotype 1 infections
# 19. new human serotype 2 infections



###############################################################
###############################################################
##### 3. Run model simulations and visualize output
###############################################################
###############################################################


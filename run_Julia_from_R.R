
###############################################################
###############################################################
##### 1. Install JuliaCall in R, load packages and environment
###############################################################
###############################################################

# install.packages("JuliaCall") 
library(ggplot2)
library(scales)
library(RColorBrewer)
library(JuliaCall)

# install_julia() #DONT RUN UNLESS YOU DON'T HAVE JULIA INSTALLED
julia <- julia_setup(JULIA_HOME = "C:/Users/ruu6/AppData/Local/Programs/Julia-1.7.2/bin")

############################################
# need "julia_..." in front of command to invoke program
# put commands in parentheses and quotation marks ("")
# if you add a semicolon at the end of the command, it will suppress the output
julia_eval("sqrt(2)")
julia_command("a = sqrt(2);")
julia_eval("a")
#julia_command will print output and return nothing, 
#julia_eval will return and print result as R object

# julia documentation: https://docs.julialang.org/en/v1/ 


############################################
# setwd("~/GitHub/sero2-model/") #set working directory to github project repo if you haven't already
# julia_install_package("Random")
# julia_install_package_if_needed("Random")
julia_library("Random,Distributions,DataFrames,LinearAlgebra,CSV,JLD2,Dates")

julia_source("1_dengue_2st.jl") #load main model function
julia_source("2_run_model_sims.jl") #load function to run main model

julia_exists("dengue_2st!") #test whether functions loaded and exist :)
julia_exists("run_model_sims!")





###############################################################
###############################################################
##### 2. Initialize parameters and populations/compartments
###############################################################
###############################################################

#####PARAMETERS
# set parameters using Wearing et al paper https://www.pnas.org/doi/full/10.1073/pnas.0602960103#:~:text=Rapid%20and%20unplanned%20urbanization%2C%20increased,contributed%20to%20dramatic%20dengue%20resurgence. 

# set seed for reproducibility
julia_command("Random.seed!(123);")

# set timeframe (in days; starts at 100 years)
tmax = julia_eval("tmax = 3650;")
julia_command("tspan = collect(0.0:(tmax+500));")

#set # simulations to run-- doing it in R
# nsims=2
nsims = julia_eval("nsims = 2;")

######## HUMAN PARAMETERS ########
julia_command("bh = 1/(60*365);") #birthrate; from paper
julia_command("m = 2;") #initial ratio of mosquitoes to humans; from literature
julia_command("c = 0.5;") #contact rate (biting rate); from paper
julia_command("p_h = 0.38;") #probability of a human being infected by an infected mosquito bite; paper 0.38
julia_command("beta_h = c*p_h*((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/365))).+ 1);") #seasonal transmission- mosquito to human
julia_command("mu_h = bh;") #human death rate; set equal to birthrate
julia_command("p_IIP = 1/5;") #progression rate out of human E state; from paper
julia_command("p_IP = 1/6;") #progression rate out of human I state; from paper
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

julia_command("par = (bh,beta_h,bm,phi_m,beta_m,mu_m,mu_mL,mu_h,p_IIP,p_IP,p_EIP,p_R,m);") #save all params into one vector to input into model function

######################################################
###########COMPARTMENTS
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

julia_command("sumu0Nm = sum(u0Sm+u0Em1+u0Im1+u0Em2+u0Im2);") #mosquito pop

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
##### 3. Run model simulations and store output
###############################################################
###############################################################

#test and run model 1 time (1 timestep)
# julia_command("x = dengue_2st!(x0,par,2)") #gives output
# x = julia_eval("dengue_2st!(x0,par,2)") #saves output as list

result = julia_eval("run_model_sims!(nsims,tmax,x0,par)") #can take a little while to run
#saved result is a 'JuliaObject' list of 34 nsims*tmax matrices 

#CASES
newcases_all_h = result[[1]] #all new human infections/cases from time 1:tmax
newcases_all_m = result[[2]] #all new mosquito infections/cases from time 1:tmax
newcases_st1_h = result[[3]] #new serotype-1 human infections/cases from time 1:tmax
newcases_st2_h = result[[4]] #new serotype-2 human infections/cases from time 1:tmax
#POPULATIONS, BIRTHS, DEATHS
hpop = result[[5]] #human population from time 1:tmax
mpop = result[[6]] #adult mosquito population from time 1:tmax
lpop = result[[7]] #larval mosquito population from time 1:tmax
hbirths = result[[8]] #human births from time 1:tmax
hdeaths = result[[9]] #human deaths from time 1:tmax

#SEROTYPE-SPECIFIC FORCES OF INFECTION
lambda_h1 = result[[10]] #mosquito-to-human serotype-1 FOI
lambda_h2 = result[[11]] #mosquito-to-human serotype-2 FOI
lambda_m1 = result[[12]] #human-to-mosquito serotype-1 FOI
lambda_m2 = result[[13]] #human-to-mosquito serotype-2 FOI

#HUMAN COMPARTMENTS
S0dt = result[[14]] #Susceptible to any serotype (0 prior infections)
E1dt = result[[15]] #Incubation compartment (IIP) after being infected by st1 for first time (only st1; no others)
E2dt = result[[16]] #Incubation compartment (IIP) after being infected by st2 for first time (only st2; no others)
I1dt = result[[17]] #Infectious compartment (IP) for st1 (only st1; no others)
I2dt = result[[18]] #Infectious compartment (IP) for st2 (only st2; no others)
R1dt = result[[19]] #Recovered and temporarily immune (cross-protection) compartment for st1 
R2dt = result[[20]] #Recovered and temporarily immune (cross-protection) compartment for st2 
S1dt = result[[21]] #Susceptible to st2 (but has been previously infected with st1)
S2dt = result[[22]] #Susceptible to st1 (but has been previously infected with st2)
E12dt = result[[23]] #Incubation compartment (IIP) after being infected by st2 (previously infected by st1)
E21dt = result[[24]] #Incubation compartment (IIP) after being infected by st1 (previously infected by st2)
I12dt = result[[25]]#Infectious compartment (IP) for st2 (previously infected by st1)
I21dt = result[[26]] #Infectious compartment (IP) for st1 (previously infected by st2)
R12dt = result[[27]] #Recovered and immune from infection with both serotypes
R21dt = result[[28]] #Recovered and immune from infection with both serotypes

#MOSQUITO COMPARTMENTS -- can only be infected once
Lmdt = result[[29]]#Development (Larval) compartment (no infections)
Smdt = result[[30]]#Susceptible to any serotype 
E1mdt = result[[31]]#Incubation compartment (EIP) after being infected by st1
E2mdt = result[[32]] #Incubation compartment (EIP) after being infected by st2
I1mdt = result[[33]]#Infectious compartment after being infected by st1
I2mdt = result[[34]] #Infectious compartment after being infected by st2



###############################################################
###############################################################
##### 4. Plot output
###############################################################
###############################################################
colors <- brewer.pal(8, name = 'Dark2')


############### NEWCASES 
#all (humans)
df <- as.data.frame(newcases_all_h)
df$median <- apply(df, 1,median)
plot(NA,NA,xlim = c(0,tmax),ylim=c(0,max(df)),ylab="New infections",xlab = "Days",bty="n")
for (i in 1:nsims){
  lines(1:tmax,df[,i],col=alpha(colors[1],0.15))
}
lines(1:tmax,df$median,col=colors[1])

#serotype 1 (humans)
df <- as.data.frame(newcases_st1_h)
df$median <- apply(df, 1,median)
plot(NA,NA,xlim = c(0,tmax),ylim=c(0,max(df)),ylab="New infections",xlab = "Days",bty="n")
for (i in 1:nsims){
  lines(1:tmax,df[,i],col=alpha(colors[2],0.15))
}
lines(1:tmax,df$median,col=colors[2])

#serotype 2 (humans)
df <- as.data.frame(newcases_st2_h)
df$median <- apply(df, 1,median)
plot(NA,NA,xlim = c(0,tmax),ylim=c(0,max(df)),ylab="New infections",xlab = "Days",bty="n")
for (i in 1:nsims){
  lines(1:tmax,df[,i],col=alpha(colors[3],0.15))
}
lines(1:tmax,df$median,col=colors[3])



############### FOI 
# serotype 1 (humans)
df <- as.data.frame(lambda_h1)
df$median <- apply(df, 1,median)
plot(NA,NA,xlim = c(0,tmax),ylim=c(0,max(df)),ylab="New infections",xlab = "Days",bty="n")
for (i in 1:nsims){
  lines(1:tmax,df[,i],col=alpha(colors[4],0.15))
}
lines(1:tmax,df$median,col=colors[4])

# serotype 2 (humans)
df <- as.data.frame(lambda_h2)
df$median <- apply(df, 1,median)
plot(NA,NA,xlim = c(0,tmax),ylim=c(0,max(df)),ylab="New infections",xlab = "Days",bty="n")
for (i in 1:nsims){
  lines(1:tmax,df[,i],col=alpha(colors[5],0.15))
}
lines(1:tmax,df$median,col=colors[5])



##### STACKED AREA PLOTS FOR COMPARTMENTS 

S.col <- brewer.pal(9, "GnBu")[7:9]
E.col <- brewer.pal(9, "Purples")[6:9]
I.col <- brewer.pal(9, "YlOrRd")[6:9]
R.col <- brewer.pal(9, "BuGn")[6:9]

source("data_stackedareaplot.R")

data <- data_stackedareaplot(S0dt,S1dt,S2dt,
                             E1dt,E2dt,E12dt,E21dt,
                             I1dt,I2dt,I12dt,I21dt,
                             R1dt,R2dt,R12dt,R21dt,
                             nsims=nsims,tmax=tmax)

#Population distribution of each compartment over timeframe
ggplot(data, aes(x=time, y=value_perc, fill=group)) + 
  geom_area() +
  theme_classic() +
  ggtitle("2-serotype model") + 
  scale_fill_manual(values=c(R.col,S.col,E.col,I.col))



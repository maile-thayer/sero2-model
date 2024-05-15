using Random,Distributions,DataFrames,LinearAlgebra,CSV,JLD2,Dates

######################################################
######################################################
# INITIALIZE MODEL
# Set seed, timeframe/timesteps, load data & variables
######################################################
######################################################

# set seed for reproducibility
Random.seed!(123)

# set timeframe
tspan = collect(0.0:50000)
tmax = 36500

nsims = 10
######################################################
# load data
pop = 3255603

######################################################
# set parameters using Wearing et al paper https://www.pnas.org/doi/full/10.1073/pnas.0602960103#:~:text=Rapid%20and%20unplanned%20urbanization%2C%20increased,contributed%20to%20dramatic%20dengue%20resurgence. 

#HUMAN
bh = 1/(60*365) #b #0.1 #per capita daily human birth rate --> 84.3 births/day islandwide, 3255642 people on island
m = 2 # initial ratio of mosquitoes to humans
c = 0.5
p_h = 0.44 #0.43 #0.4
# p_h = 0.38 # paper
βh = c*p_h*((0.3 *(cos.(((2*pi*tspan).+ 5.295594)/365))).+ 1) 
μh = bh

#MOSQUITO
bm = 1/14 #repeat([0.2],n) #mosquito laying rate #0.06
ψm = 1 #ones(n) #emergence rate
p_m = p_h
βm = c*p_m
μm = bm # adult mosquito death rate
μm_L = 0.0  #death rate for larval stage

# probabilities for leaving gamma-distributed sub-compartments
p_EIP = 1/10
p_IIP = 1/5
p_IP = 1/6
p_R = 1/365 #or 1/120 

par = [bh,βh,bm,ψm,βm,μm,μm_L,μh,p_IIP,p_IP,p_EIP,p_R]; #save all single parameters in one vector

R0 = (m*median(βh)*βm*p_EIP)/(μm*(p_IP+bh)*(p_EIP+μm))


######################################################
# INITIALIZE POPULATIONS
u0Sh_0 = round(0.423*pop) 
# u0Ih_1 = ones(n,a,15)*(1)
u0Eh_1 = 0 #all Es are empty to start
u0Ih_1 = 25 #start with 100 infectious in I1
u0Rh_1 = round((2/70)*0.227*pop)

u0Eh_2 = 0
u0Ih_2 = 25
u0Rh_2 = round((2/70)*0.227*pop)

u0Sh_1 = round((68/70)*0.227*pop)
u0Eh_12 = 0
u0Ih_12 = 25
u0Rh_12 = round(0.061*pop)

u0Sh_2 = round((68/70)*0.227*pop)
u0Eh_21 = 0
u0Ih_21 = 25
u0Rh_21 = round(0.061*pop) 


u0Lm = 0
u0Sm = pop*m
u0Em1 = 0
u0Im1 = 0
u0Em2 = 0
u0Im2 = 0

sumu0Nm = sum(u0Sm+u0Em1+u0Im1+u0Em2+u0Im2) 

# Vector of all compartments
u0all = ([u0Sh_0, u0Eh_1, u0Eh_2, u0Ih_1, u0Ih_2, u0Rh_1, u0Rh_2, u0Sh_1, u0Sh_2, u0Eh_12, u0Eh_21, u0Ih_12, u0Ih_21, u0Rh_12, u0Rh_21, u0Lm, u0Sm, u0Em1, u0Em2, u0Im1, u0Im2])

intv = "nointv"
# list of all initial items of interest AT TIME 0, mimicking model input and output
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
# X 13. aging (by age) 
# X 14. deaths by age
# X 15. importations
# 16. new secondary human cases
# 17. total larval population
# 18. new human serotype 1 infections
# 19. new human serotype 2 infections
x0 = [u0all,0,0,pop,sumu0Nm,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
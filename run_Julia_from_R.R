

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

julia_source("1_Functions_sero2.jl") # load main model function
julia_source("2_Initialization_sero2.jl")
julia_source("3_Runmodel_sero2.jl")

julia_exists("dengue_sde_2st!")

# using Random,Distributions,DataFrames,LinearAlgebra,CSV,JLD2,Dates
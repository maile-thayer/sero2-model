

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




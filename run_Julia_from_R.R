

install.packages("JuliaCall")

library(JuliaCall)
install_julia()
julia <- julia_setup()

############################################
# need "julia_..." in front of command to invoke program
# put commands in parenthases ""
#test out sqrt(2)
julia_eval("sqrt(2)")
#assign
julia_command("a = sqrt(2);")
julia_eval("a")

############################################
setwd("~/GitHub/DE_model/")


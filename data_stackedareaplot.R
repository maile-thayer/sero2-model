data_stackedareaplot <- function(S0dt,S1dt,S2dt,
                                 E1dt,E2dt,E12dt,E21dt,
                                 I1dt,I2dt,I12dt,I21dt,
                                 R1dt,R2dt,R12dt,R21dt,nsims=10,tmax=36500){
  #add row median column
  S0dt = cbind(S0dt,apply(E1dt, 1,median)) 
  E1dt = cbind(E1dt,apply(E1dt, 1,median)) 
  E2dt = cbind(E2dt,apply(E2dt, 1,median)) 
  I1dt = cbind(I1dt,apply(I1dt, 1,median)) 
  I2dt = cbind(I2dt,apply(I2dt, 1,median)) 
  R1dt = cbind(R1dt,apply(R1dt, 1,median)) 
  R2dt = cbind(R2dt,apply(R2dt, 1,median)) 
  S1dt = cbind(S1dt,apply(S1dt, 1,median)) 
  S2dt = cbind(S2dt,apply(S2dt, 1,median)) 
  E12dt = cbind(E12dt,apply(E12dt, 1,median)) 
  E21dt = cbind(E21dt,apply(E21dt, 1,median)) 
  I12dt = cbind(I12dt,apply(I12dt, 1,median)) 
  I21dt = cbind(I21dt,apply(I21dt, 1,median)) 
  R12dt = cbind(R12dt,apply(R12dt, 1,median)) 
  R21dt = cbind(R21dt,apply(R21dt, 1,median)) 
  
  foo <- (as.data.frame((rbind(S0dt[,(nsims+1)],S1dt[,(nsims+1)],S2dt[,(nsims+1)],
                               E1dt[,(nsims+1)],E2dt[,(nsims+1)],E12dt[,(nsims+1)],E21dt[,(nsims+1)],
                               I1dt[,(nsims+1)],I2dt[,(nsims+1)],I12dt[,(nsims+1)],I21dt[,(nsims+1)],
                               R1dt[,(nsims+1)],R2dt[,(nsims+1)],R12dt[,(nsims+1)],R21dt[,(nsims+1)]))))
  
  N <- colSums(foo)
  
  
  time <- as.numeric(rep(seq(1,tmax),times=15))
  # value_count <- as.vector(t(as.matrix(foo)))
  value_perc <- c(S0dt[,(nsims+1)]/N,
                  S1dt[,(nsims+1)]/N,
                  S2dt[,(nsims+1)]/N,
                  E1dt[,(nsims+1)]/N,
                  E2dt[,(nsims+1)]/N,
                  E12dt[,(nsims+1)]/N,
                  E21dt[,(nsims+1)]/N,
                  I1dt[,(nsims+1)]/N,
                  I2dt[,(nsims+1)]/N,
                  I12dt[,(nsims+1)]/N,
                  I21dt[,(nsims+1)]/N,
                  R1dt[,(nsims+1)]/N,
                  R2dt[,(nsims+1)]/N,
                  R12dt[,(nsims+1)]/N,
                  R21dt[,(nsims+1)]/N)
  
  group <- rep(c("S0h","S1h","S2h",
                 "E1h","E2h","E12h","E21h",
                 "I1h","I2h","I12h","I21h",
                 "R1h","R2h","R12h","R21h"),each=tmax)        
  data <- data.frame(time, value_perc, group)
  data$group <- factor(data$group, levels = c("R21h","R12h","R2h","R1h",
                                              "S2h","S1h","S0h",
                                              "E21h","E12h","E2h","E1h",
                                              "I12h","I21h","I2h","I1h"
  ))
  return(data)
  
}



## Simulation of the SIR Model ====

# load packages
library(cowplot)
library(deSolve)
library(ggplot2)
library(graphics)
library(reshape2)
library(tidyverse)

# clear environment
rm(list=ls()) 


## Initialize ====

# time 
dt = 0.05
tvals = seq(0,150,by=dt)

# parameters
pars = c(beta=0.00009, gamma=0.05)

# initial conditions
N = 10000
inits = c(S=N-1, I=1, R=0, Cases=0)


## Model ====

SIR = function(t, inits, pars){
  with(as.list(c(inits, pars)), {
    dS = -beta*S*I
    dI = beta*S*I-gamma*I
    dR = gamma*I
    dCases = beta*S*I
    list(c(dS, dI, dR, dCases))
  })
}

# solve model
sol = as.data.frame(ode(y=inits, times=tvals, func=SIR, parms=pars))


## Total Cases ====

sol$Cases[length(sol$Cases)]

# peak
max(sol$I)
which.max(sol$I)*dt


## Plot ====

sol %>% gather(key, individuals, S, I, R) %>%
  ggplot(aes(x=time, y=individuals, color=key)) +
  geom_line(size=1.2) +
  ggtitle("SIR model") + 
  xlab("days") + ylab("individuals") +
  scale_color_manual(values=c("#2623D6","#D62323","#468220"),
                     name="key", 
                     breaks=c("S","I","R"),
                     labels=c("S","I","R")) +
  theme(text=element_text(size=15),
        legend.position=c(1,1),
        legend.justification=c("right","top"),
        legend.margin=margin(5,5,5,5))

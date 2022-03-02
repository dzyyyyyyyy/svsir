library(deSolve)
library(ggplot2)

# initial state explanations:
#   -VjN = proportion of population who have received vaccine N
#   -S = proportion of population susceptible to infection, non-vaccinated
#   -I = proportion of population currently infected with disease
#   -R = proportion of population recovered

Vj1 = 0
Vj2 = 0
Vj3 = 0
I = 0.1
R = 0
S = 1- Vj1-Vj2-Vj3-I-R
N = Vj1 + Vj2 + Vj3 + S + I + R 

initial_state_values = c(Vj1 = Vj1, Vj2 = Vj2, Vj3 = Vj3, S = S, 
                         I = I, R= R)



# parameter explanations:
#   - beta: infection rate / day
#   - gamma: recovery rate / day
#   - u: birth/death rate
#   - elN: leakage from vaccine N
#   - eaN: primary vaccine failure rate for vaccine N (after you get the vaccine, how many people lose immunity)
#   - ajN: waning rate for vaccine N (rate of people which don't receive any protection from vaccine). here, set to zero
#   - rjN: rate of vaccination for vaccine N 
#   - doses: number of doses available for initial use
#   - population: population in the country


#these particular parameters are outside the vector  bc i think i had some trouble w/ exponential notation 
#within the named vector, so its easier to just tell R it's a number here then reference that in the parameters vector

rj1 = 0.02
rj2 = 0.02
rj3 = 0

doses = 5e5
population = 1e6

parameters = c(beta = 0.11808, gamma = 0.09449, u =0.00023, 
               el1 = 0.8, ea1 = 3.77e-4, aj1 = 0, rj1 = rj1,
               el2 = 0.13, ea2 = 3.77e-4, aj2 = 0, rj2 = rj2,
               el3 = 0.002, ea3 =0.00268, aj3 = 0, rj3 = rj3,
               doses = doses , population = population)

time = seq(from=1,to=200,by=1)

vsir_model <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    N = Vj1 + Vj2 + Vj3 + S + I + R
    
    rj2t = rj2 # save the rj2 variable for later
    rj3t = rj3
    
    if (time < 30){
      rj2 = 0 # we set back rj2 to zero until 30 days have passed (time you have to wait from 1st shot to 2nd for COVID vaccine)
    } else {
      rj2 = rj2t # set rj2 to rj2t value
    }
    
    
    # this expression isn't great bc it will assume that people who are initially dosed
    # in your parameters got their doses from the supply you want to start applying
    # but for our question (what's the best thing to do w/ 0 people vaxxed), this should work fine
    # also it'll go over the amount of doses you actually have for one iteration. not great
    # this also makes the ode really slow and increases maxsteps somehow. man this code sucks
    
    doses_left = doses - (Vj1 + 2*Vj2 + 3*Vj3)*population 
    doses = doses_left
    
    if (doses_left <= 1){
      rj1 = 0
      rj2 = 0
      rj3 = 0
      rj2t = 0 
      rj3t = 0
    }

    # change in vaccinated compartment
    dV1 = (1-ea1)*(rj1)*S - (aj1 + u)*Vj1 - el1*beta*I*Vj1 - (1-ea2)*rj2*Vj1
    dV2 = (1-ea2)*(rj2)*Vj1 - (aj2 + u)*Vj2 - el2*beta*I*Vj2 - (1-ea3)*rj3*Vj2
    dV3 = (1-ea3)*(rj3)*Vj2 - (aj3 + u)*Vj3 - el3*beta*I*Vj3
    
    
    
    # change in susceptible compartment
    dS = u + 
      (aj1*Vj1-(1-ea1)*rj1*S) +
      (aj2*Vj2) + # took out other term here -- people are not going from susceptible to Vj2
      (aj3*Vj3) - # took out other term here -- people are not going from susceptible to Vj3
      beta*S*I - 
      u * S
    
    # change in infected compartment
    dI = beta*I*
      ((el1*Vj1)+(el2*Vj2)+(el3*Vj3))+
      beta*S*I -
      (gamma+u)*I
    
    # change in recovered compartment
    dR = gamma*I - u*R
    
    return(list(c(dV1, dV2, dV3, dS,dI,dR)))
  }
  )
}

output<-as.data.frame(ode(y=initial_state_values,func = vsir_model,parms=parameters,times = time))
output$N <- output$Vj1 + output$Vj2 + output$Vj3 + output$S + output$I + output$R


out_long = reshape2::melt(output, id ='time')

ggplot(data = out_long,
       aes(x=time, y=value, colour = variable, group = variable)) +
  geom_line() +xlab('Time (days') + ylab('Propotion of the population') + scale_color_discrete(name="State")


output$doses_used <- (output$Vj1 + 2*output$Vj2 + 3*output$Vj3)*population 


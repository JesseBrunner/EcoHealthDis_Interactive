library(tidyverse)
library(simecol)


# The model
SEIR <- odeModel(
  main = function (time, init, parms, ...) {
    # (1) unpack variables
    S <- init["S"]
    E <- init["E"]
    I <- init["I"]
    R <- init["R"]
    C <- init["C"]

    N <- S+E+I+R

    with(as.list(parms),  {
      # Model from Chowell et al. 2004 The basic reproductive number of Ebola and the effects of public health measures: the cases of Congo and Uganda. Journal of Theoretical Biology 229:119-126.

      dS <- -beta*S*I/N
      dE <-  beta*S*I/N  -k*E
      dI <-               k*E -gamma*I
      dR <-                    gamma*I
      dC <-               k*E
      list(c(dS, dE, dI, dR, dC))
    })
  },
  solver = "rk4" # the function that does the solving
)

# Congo data set
Congo <- read_csv("Ebola_Congo_1995.csv"
)
Congo <- Congo %>%
  mutate(Cumulative = cumsum(Cases)) %>%
  filter(Day >=7)

parms_congo <- c(beta=0.33, k=1/5.3, gamma=1/5.61)
inits_congo <- c(S=5364500, E=1, I=0, R=0, C=1)
times_congo = c(from=-7, to=128, by=1) #by=step size


# Uganda data set
Uganda <- read_csv("Ebola_Uganda_2000.csv")
Uganda <- Uganda %>%
  mutate(Cumulative = cumsum(Cases))
parms_uganda <- c(beta=0.4, k=1/3.35, gamma=1/3.5)
inits_uganda <- c(S=1867200, E=0, I=3, R=4, C=7)
times_uganda = c(from=0, to=133, by=1) #by=step size

# choose data set
df <- Congo
pars <- parms_congo
inits <- inits_congo
timespan <- times_congo

df <- Uganda
pars <- parms_uganda
inits <- inits_uganda
timespan <- times_uganda

# set parameters
parms(SEIR) <- pars
init(SEIR) <- inits
times(SEIR) <- timespan

SEIR <- sim(SEIR)

out(SEIR) %>% # get the output from the model
  filter(time >=0) %>%
  ggplot(., aes(x=time, y=C)) + # construct the plot
  geom_line() +
  geom_point(data=df, aes(Day, Cumulative)) + scale_y_log10()

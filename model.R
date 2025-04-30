library(deSolve)
library(ggplot2)
library(tidyr)


eggs<-seq(0,1000,by=10)
control<-0#seq(0,.5,by=.1)#seq(0,1,by=0.5)
outputlist<-list()
index<-1
Maxbiomass<-50
HarvestIndex<-0.2
for(l in eggs){
  for(c in control){
# Parameters from the table
params <- c(        #Initial Mass
  r = 0.1,        # Growth rate of maize
  K_M = Maxbiomass,        # Max biomass of maize
  beta = 5e-8,     # Plant attack rate by caterpillars
  b = 125,         # Eggs laid per day per female
  W = 0.5,         # Proportion of females
  K_E = 1e8,       # Carrying capacity of eggs
  K_L = 1e6,       # Carrying capacity of larvae
  alpha_E = 1/3,   # Rate out of egg stage
  mu_E = 0.01,     # Mortality of eggs
  alpha_L = 1/14,  # Rate out of larvae
  mu_L = 0.01,     # Mortality of larvae
  alpha_P = 1/9,   # Rate out of pupae
  mu_P = 0.01,     # Mortality of pupae
  mu_A = 0.01,     # Mortality of adults
  alpha_p = 1/9,   # Same as alpha_P
  e = 0.2,          # Leaf impact factor
  u_E = c,          # Egg control factor
  u_L = c,         # Larvae control factor
  u_P = c,          # Pupae control factor
  u_A = c          # Adult control factor
)

# Derived parameter
params["theta"] <- params["e"] * params["beta"]

# Initial conditions: M, E, L, P, A
state <- c(M = 15, E = l, L = 0, P = 0, A = 0)

# Time
times <- seq(0, 120, by = 1)

# ODE system
faw_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dM <- (r * M * (1 - (M / K_M))) - (theta * L * M)
    dE <- (b * (1 - (E / K_E)) * W * A) - ((alpha_E + mu_E + u_E) * E)
    dL <- (alpha_E * (1 - (L / K_L)) * E) + (theta * L * M) - ((alpha_L + mu_L + u_L) * L)
    dP <- (alpha_L * L) - ((mu_P + alpha_p + u_P) * P)
    dA <- (alpha_p * P) - ((mu_A + u_A) * A)
    
    list(c(dM, dE, dL, dP, dA))
  })
}

# Solve the system
out <- ode(y = state, times = times, func = faw_model, parms = params)

# Convert to data frame
outputlist[[index]] <- data.frame(Biomass = max(out[,2]), Control_Proportion = c,
                                  Eggs = l)
index<-index+1
}
}
outdf<-do.call(rbind,outputlist)
outdf$Yieldestimate<-outdf$Biomass*HarvestIndex
outdf$Yieldloss<-100-((outdf$Yieldestimate/(HarvestIndex*Maxbiomass) *100))
p <- ggplot(outdf, aes(x = Eggs, y = Yieldloss)) +
  geom_line(size = 1) +
  labs(
    x = "Initial Egg Density",
    y = "Yield Loss (%)"
   ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), # for export
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.position = "right",
    legend.key = element_blank()
  )


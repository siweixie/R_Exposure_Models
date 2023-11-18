# Define model 200
model_200 <- function(G, Q, beta, gamma) {
  C_F <- gamma * G / Q # Far field
  C_N <- C_F + ((gamma * G) / beta) # Near field
  return(list(C_F = C_F, C_N = C_N))
}

# Define model 201
# Define Model 105
model_201 <- function(G, Q, Q_L, epsilon_L, V, t_g, T) {
  # Create a time vector
  time_vector <- seq(0, T, by = 1)
  
  C_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return(((G * (1 - epsilon_L)) / (Q + Q_L)) * (1 - exp((-(Q + Q_L) * t) / V)))
    } else {
      return(NA)
    }
  })
  
  C0 <- (G * (1 - epsilon_L)) / (Q + Q_L) * (1 - exp((-(Q + Q_L) * t_g) / V))
  
  C_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(C0 * exp((-(Q + Q_L) * (t - t_g)) / V))
    } else {
      return(NA)
    }
  })
  
  return(list(time = time_vector, concentration_rise = C_rise, concentration_decay = C_decay))
}


# Define the model parameters with units
td <- 60    # minutes
tg <- 15    # minutes
G <- 100    # mg/min
Q <- 20     # m^3/min
QL <- 5     # m^3/min
epsilon_L <- 0.5  # Local control efficiency
epsilon_LF <- 0.75 # Filtration efficiency
Q_R <- 5    # m^3/min
V <- 100    # m^3
V_N <- 8    # Near-field volume in m^3
beta <- 5   # Near-field ventilation rate in m^3/min
gamma <- 0.25


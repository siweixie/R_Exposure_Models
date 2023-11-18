# Define Model 104
model_104 <- function(G, Q, Q_L, epsilon_L, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L)) / (Q + Q_L)
  C_L_avg <- C_avg + (epsilon_L * gamma * G) / Q_L
  return(list(C_avg = C_avg, C_L_avg = C_L_avg))
}


# Define Model 105
model_105 <- function(G, Q, Q_L, epsilon_L, V, t_g, T) {
  time_vector <- seq(0, T, by = 1)
  
  C_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return((G * (1 - epsilon_L)) / (Q + Q_L) * (1 - exp((-(Q + Q_L) * t) / V)))
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

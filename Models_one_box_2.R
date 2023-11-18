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


# Define Model 106
model_106 <- function(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L)) / ((Q + epsilon_RF * Q_R) + Q_L)
  C_L_avg <- C_avg + (epsilon_L * gamma * G) / Q_L
  C_RF <- C_avg * (1 - epsilon_RF)
  return(list(C_avg = C_avg, C_L_avg = C_L_avg, C_RF = C_RF))
}


model_107 <- function(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, V, t_g, T) {
  Q_instead <- Q + epsilon_RF * Q_R
  return(model_105(G, Q_instead, Q_L, epsilon_L, V, t_g, T))
}


results_104 <- model_104(G, Q, Q_L, epsilon_L, gamma)
results_105 <- model_105(G, Q, Q_L, epsilon_L, V, t_g, T)
results_106 <- model_106(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, gamma)
results_107 <- model_107(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, V, t_g, T)

results_105_df <- data.frame(Time = results_105$time, 
                             Concentration = c(results_105$concentration_rise, results_105$concentration_decay))
results_105_df <- results_105_df[!is.na(results_105_df$Concentration), ]

results_107_df <- data.frame(Time = results_107$time, 
                             Concentration = c(results_107$concentration_rise, results_107$concentration_decay))
results_107_df <- results_107_df[!is.na(results_107_df$Concentration), ]



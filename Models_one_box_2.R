# Define Model 104
model_104 <- function(G, Q, Q_L, epsilon_L, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L)) / (Q + Q_L)
  C_L_avg <- C_avg + (epsilon_L * gamma * G) / Q_L
  return(list(C_avg = C_avg, C_L_avg = C_L_avg))
}


# Define Model 105
model_105 <- function(G, Q, Q_L, epsilon_L, V, t_g, T) {
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


# Define Model 106
model_106 <- function(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L)) / ((Q + epsilon_RF * Q_R) + Q_L)
  C_L_avg <- C_avg + (epsilon_L * gamma * G) / Q_L
  C_RF <- C_avg * (1 - epsilon_RF)
  return(list(C_avg = C_avg, C_L_avg = C_L_avg, C_RF = C_RF))
}


# Define Model 107
model_107 <- function(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, V, t_g, T) {
  Q_instead <- Q + epsilon_RF * Q_R
  return(model_105(G, Q_instead, Q_L, epsilon_L, V, t_g, T))
}


T <- 60   # Total time (minutes)
t_g <- 15 # Time of generation (minutes)
G <- 100  # mg/min
Q <- 20   # m^3/min
Q_L <- 5  # m^3/min
epsilon_L <- 0.5  # Efficiency of local exhaust
epsilon_L_F <- 0.75  # Efficiency of local exhaust filtration
Q_R <- 5  # m^3/min
epsilon_RF <- 0.9  # Efficiency of recirculation filtration
V <- 100  # m^3
gamma <- 0.25 


results_104 <- model_104(G, Q, Q_L, epsilon_L, gamma)
results_105 <- model_105(G, Q, Q_L, epsilon_L, V, t_g, T)
results_106 <- model_106(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, gamma)
results_107 <- model_107(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, V, t_g, T)

results_105_df <- data.frame(Time = results_105$time, 
                             Concentration = c(results_105$concentration_rise, 
                                               results_105$concentration_decay,
                                              Model = "Model 105"))
results_105_df <- results_105_df[!is.na(results_105_df$Concentration), ]

results_107_df <- data.frame(Time = results_107$time, 
                             Concentration = c(results_107$concentration_rise,
                                               results_107$concentration_decay,
                                              Model = "Model 107"))
results_107_df <- results_107_df[!is.na(results_107_df$Concentration), ]


# Comparison of model 105 and model 107:

combined_data <- rbind(results_105_df, results_107_df)
ggplot(combined_data, aes(x = Time, y = Concentration, color = Model)) +
  geom_line() +
  labs(title = "Model 101 vs Model 103",
       x = "Time (minutes)", 
       y = "Concentration (mg/m^3)") +
  theme_minimal()


# Further check the significance of difference
mse <- mean((results_105_df$Concentration - results_107_df$Concentration)^2)
r2 <- cor(results_105_df$Concentration, results_107_df$Concentration)^2
mse
r2

# Normality test
shapiro.test(results_105_df$Concentration)
shapiro.test(results_107_df$Concentration)

# If it is not normal distribution (p<0.05), using wilcox.test
wilcox.test(results_105_df$Concentration, results_107_df$Concentration)

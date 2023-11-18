# Define Model 108
model_108 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)
  C_LE <- C_avg + (epsilon_L * gamma * G) / Q_L
  C_LF <- C_LE * (1 - epsilon_L_F)
  return(list(C_avg = C_avg, C_LE = C_LE, C_LF = C_LF))
}


# Define Model 109
model_109 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, V, t_g, T) {
  # Create a time vector
  time_vector <- seq(0, T, by = 1)
  
  C_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return(((G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)) * (1 - exp((-(Q + epsilon_L_F * Q_L) * t) / V)))
    } else {
      return(NA)
    }
  })
  
  C0 <- ((G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)) * (1 - exp((-(Q + epsilon_L_F * Q_L) * t_g) / V))
  
  C_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(C0 * exp((-(Q + epsilon_L_F * Q_L) * (t - t_g)) / V))
    } else {
      return(NA)
    }
  })
  
  return(list(time = time_vector, concentration_rise = C_rise, concentration_decay = C_decay))
}


# Define Model 110
model_110 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L * epsilon_L_F)) / ((Q + epsilon_RF * Q_R) + epsilon_L_F * Q_L)
  C_LE <- C_avg + (epsilon_L * gamma * G) / Q_L
  C_LF <- C_LE * (1 - epsilon_L_F)
  C_RF <- C_avg * (1 - epsilon_RF)
  return(list(C_avg = C_avg, C_LE = C_LE, C_LF = C_LF, C_RF = C_RF))
}


# Define Model 111
model_111 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, V, t_g, T) {
  Q_instead <- Q + epsilon_RF * Q_R
  return(model_109(G, Q_instead, Q_L, epsilon_L, epsilon_L_F, V, t_g, T))
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


results_108 <- model_108(G, Q, Q_L, epsilon_L, epsilon_L_F, gamma)
results_109 <- model_109(G, Q, Q_L, epsilon_L, epsilon_L_F, V, t_g, T)
results_110 <- model_110(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, gamma)
results_111 <- model_111(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, V, t_g, T)

results_109_df <- data.frame(Time = results_109$time, 
                             Concentration = c(results_109$concentration_rise, 
                                               results_109$concentration_decay),
                                              Model = "Model 109")
results_109_df <- results_109_df[!is.na(results_109_df$Concentration), ]

results_111_df <- data.frame(Time = results_111$time, 
                             Concentration = c(results_111$concentration_rise,
                                               results_111$concentration_decay),
                                              Model = "Model 111")
results_111_df <- results_111_df[!is.na(results_111_df$Concentration), ]


# Comparison of model 109 and model 111:

combined_data <- rbind(results_109_df, results_111_df)
ggplot(combined_data, aes(x = Time, y = Concentration, color = Model)) +
  geom_line() +
  labs(title = "Model 109 vs Model 111",
       x = "Time (minutes)", 
       y = "Concentration (mg/m^3)") +
  theme_minimal()


# Further check the significance of difference
mse <- mean((results_109_df$Concentration - results_111_df$Concentration)^2)
r2 <- cor(results_109_df$Concentration, results_111_df$Concentration)^2
mse
r2

# Normality test
shapiro.test(results_109_df$Concentration)
shapiro.test(results_111_df$Concentration)

# If it is not normal distribution (p<0.05), using wilcox.test
wilcox.test(results_109_df$Concentration, results_111_df$Concentration, paired = TRUE)


# Overall comparison
combined_data <- rbind(results_101_df, results_103_df, results_105_df, results_107_df, results_109_df, results_111_df)
ggplot(combined_data, aes(x = Time, y = Concentration, color = Model)) +
    geom_line() +
    labs(title = "Model difference",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal()

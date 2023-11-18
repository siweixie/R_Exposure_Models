# Define Model 100 with the gamma factor
model_100 <- function(G, Q, gamma) {
  C_steady <- gamma * G / Q
  return(C_steady)
}

# Define Model 101
model_101 <- function(G, Q, V, t_g, T) {
  # Create a time vector
  time_vector <- seq(0, T, by = 1)
  
# Calculate the phase of concentration increase
  C_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return((G / Q) * (1 - exp((-Q * t) / V )))
    } else {
      return(NA)
    }
  })          
  
  # Concentration at the end of the emission phase as the initial concentration for the decay phase
  C0 <- (G / Q) * (1 - exp(-Q * t_g / V ))
  
  # Calculate the concentration decay phase
  C_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(C0 * exp(-Q * (t - t_g) / V ))
    } else {
      return(NA)
    }
  })
  
  return(list(time = time_vector, concentration_rise = C_rise, concentration_decay = C_decay))
}
             

# Calculation of steady-state concentrations for model 100
C_steady_100 <- model_100(G, Q, gamma)

# Calculation of steady-state concentrations for model 101
results_101 <- model_101(G, Q, V, t_g, T)


# Define Model 102
model_102 <- function(G, Q, Q_R, epsilon_RF, gamma) {
  C_bar <- (gamma * G) / (Q + epsilon_RF * Q_R)
  C_RF <- C_bar * (1 - epsilon_RF)
  return(list(C_bar = C_bar, C_RF = C_RF))
}


# Define Model 103
model_103 <- function(G, Q, Q_R, epsilon_RF, V, t_g, T) {
  Q_instead <- Q + epsilon_RF * Q_R
  results_3 <- model_101(G, Q_instead, V, t_g, T)
  return(results_3)
}

T <- 60   # Total time (minutes)
t_g <- 15 # Time of generation (minutes)
G <- 100  # mg/min
Q <- 20   # m^3/min
epsilon_L <- 0.5  # Efficiency of local exhaust
epsilon_L_F <- 0.75  # Efficiency of local exhaust filtration
Q_R <- 5  # m^3/min
epsilon_RF <- 0.9  # Efficiency of recirculation filtration
V <- 100  # m^3
gamma <- 0.25  

results_102 <- model_102(G, Q, Q_R, epsilon_RF, gamma)
results_103 <- model_103(G, Q, Q_R, epsilon_RF, V, t_g, T)

#Visualization
library(ggplot2)
C_bar_1 <- results_102$C_bar

results_103_df <- data.frame(Time = results_103$time, 
                             Concentration = c(results_103$concentration_rise, results_103$concentration_decay))
results_103_df <- results_103_df[!is.na(results_103_df$Concentration), ]

ggplot() + 
    geom_line(data = results_103_df, aes(x = Time, y = Concentration), colour = "blue") +
    geom_hline(yintercept = C_bar_1, linetype = "dashed", color = "red") +
    labs(title = "Model 102 vs Model 103", 
         x = "Time (minutes)", y = "Concentration (mg/m^3)") +
    theme_minimal()


# Comparison of model 100 and model 102:
comparison_0_2 <- data.frame(Model = c("Model 100", "Model 102"),
                            Concentration = c(C_steady_100, results_102$C_bar))

ggplot(comparison_0_2, aes(x = Model, y = Concentration, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "Model 100 vs Model 102",
       x = "Model",
       y = "Steady State Concentration (mg/m^3)") +
  theme_minimal()


# Comparison of model 101 and model 103:
results_101_df <- data.frame(Time = results_101$time, 
                             Concentration = c(results_101$concentration_rise, results_101$concentration_decay),
                            Model = "Model 101")

results_101_df <- results_101_df[!is.na(results_101_df$Concentration), ]

results_103_df <- data.frame(Time = results_103$time, 
                             Concentration = c(results_103$concentration_rise, results_103$concentration_decay),
                            Model = "Model 103")
results_103_df <- results_103_df[!is.na(results_103_df$Concentration), ]

combined_data <- rbind(results_101_df, results_103_df)
ggplot(combined_data, aes(x = Time, y = Concentration, color = Model)) +
  geom_line() +
  labs(title = "Model 101 vs Model 103",
       x = "Time (minutes)", 
       y = "Concentration (mg/m^3)") +
  theme_minimal()


# Further check the significance of difference
mse <- mean((results_101_df$Concentration - results_103_df$Concentration)^2)
r2 <- cor(results_101_df$Concentration, results_103_df$Concentration)^2
mse
r2

# Normality test
shapiro.test(results_101_df$Concentration)
shapiro.test(results_103_df$Concentration)

# If it is not normal distribution (p<0.05), using wilcox.test
wilcox.test(results_101_df$Concentration, results_103_df$Concentration)



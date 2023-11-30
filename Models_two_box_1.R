library(ggplot2)

# Define model 200
model_200 <- function(G, Q, beta, gamma) {
  C_F <- gamma * G / Q # Far field
  C_N <- C_F + ((gamma * G) / beta) # Near field
  return(list(C_F = C_F, C_N = C_N))
}

# Define model 201

# Define the model parameters with units
T <- 60    # minutes
t_g <- 15    # minutes
G <- 100    # mg/min
Q <- 20     # m^3/min
Q_L <- 5     # m^3/min
epsilon_L <- 0.5  # Local control efficiency
epsilon_LF <- 0.75 # Filtration efficiency
epsilon_RF <- 0.9
Q_R <- 5    # m^3/min
V <- 100    # m^3
V_N <- 8    # Near-field volume in m^3
beta <- 5   # Near-field ventilation rate in m^3/min
gamma <- 0.25

# Derived parameters
V_F <- V - V_N  # Far-field volume in m^3

# Time points
time_points <- 0:T

# Lambda calculations
lambda_1 <- 0.5 * (-(beta * V_F + V_N * (beta + Q)) / (V_N * V_F) +
                     sqrt(((beta * V_F + V_N * (beta + Q)) / (V_N * V_F))^2 -
                            4 * (beta * Q) / (V_N * V_F)))
lambda_2 <- 0.5 * (-(beta * V_F + V_N * (beta + Q)) / (V_N * V_F) -
                     sqrt(((beta * V_F + V_N * (beta + Q)) / (V_N * V_F))^2 -
                            4 * (beta * Q) / (V_N * V_F)))

# Initial concentrations for the decay phase
C_F0 <- G / Q  # Steady state far field concentration
C_N0 <- C_F0 + G / beta  # Steady state near field concentration

# Concentration calculations
C_F_t <- numeric(length(time_points))
C_N_t <- numeric(length(time_points))

for (t in time_points) {
  if (t <= t_g) {
    C_F_t[t + 1] <- C_F0 + G * ((lambda_1 * V_N + beta) / beta) * 
      ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
      exp(lambda_1 * t) - G * ((lambda_2 * V_N + beta) / beta) * 
      ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
      exp(lambda_2 * t)
    C_N_t[t + 1] <- C_N0 + G * ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
      exp(lambda_1 * t) - G * 
      ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
      exp(lambda_2 * t)
  } else {
    C_F_t[t + 1] <- C_F_t[t_g + 1] * exp(lambda_1 * (t - t_g))
    C_N_t[t + 1] <- C_N_t[t_g + 1] * exp(lambda_1 * (t - t_g))
  }
}


results_201 <- data.frame(Time = time_points, C_F_t = C_F_t, C_N_t = C_N_t)


# Model 202
model_202 <- function(G, Q, Q_R, beta, gamma, epsilon_RF) {
  C_F <- gamma * G / (Q + Q_R * epsilon_RF) # Steady state far field concentration
  C_N <- C_F + (gamma * G) / beta           # Steady state near field concentration
  C_R_F <- C_F * (1 - epsilon_RF)           # Concentration on the return side after filtration
  
  return(list(C_F = C_F, C_N = C_N, C_R_F = C_R_F))
}

# Run Model 202
results_202 <- model_202(G, Q, Q_R, beta, gamma, epsilon_LF)



# Model 203
model_203 <- function(G, Q, Q_R, beta, gamma, epsilon_RF, epsilon_L, t_g, T, V, V_N) {
  # Adjusted ventilation rate including filtered recirculated air
  Q_prime <- Q + epsilon_RF * Q_R
  
  # Derived parameters
  V_F <- V - V_N  # Far-field volume in m^3
  
  # Time points
  time_points <- 0:T
  
  # Lambda calculations (same as in Model 201 but with Q replaced by Q_prime)
  lambda_1 <- 0.5 * (-(beta * V_F + V_N * (beta + Q_prime)) / (V_N * V_F) +
                       sqrt(((beta * V_F + V_N * (beta + Q_prime)) / (V_N * V_F))^2 -
                              4 * (beta * Q_prime) / (V_N * V_F)))
  lambda_2 <- 0.5 * (-(beta * V_F + V_N * (beta + Q_prime)) / (V_N * V_F) -
                       sqrt(((beta * V_F + V_N * (beta + Q_prime)) / (V_N * V_F))^2 -
                              4 * (beta * Q_prime) / (V_N * V_F)))
  
  
  # Initial concentrations for the decay phase with Q replaced by Q_prime
  C_F0 <- G / Q_prime  # Steady state far field concentration
  C_N0 <- C_F0 + G / beta  # Steady state near field concentration
  
  # Concentration calculations
  C_F_t <- numeric(length(time_points))
  C_N_t <- numeric(length(time_points))
  
  for (t in time_points) {
    if (t <= t_g) {
      # Concentration rise equations with Q replaced by Q_prime
      C_F_t[t + 1] <- C_F0 + G * ((lambda_1 * V_N + beta) / beta) * 
        ((beta * Q_prime + lambda_2 * V_N * (beta + Q_prime)) / (beta * Q_prime * V_N * (lambda_1 - lambda_2))) * 
        exp(lambda_1 * t) - G * ((lambda_2 * V_N + beta) / beta) * 
        ((beta * Q_prime + lambda_1 * V_N * (beta + Q_prime)) / (beta * Q_prime * V_N * (lambda_1 - lambda_2))) * 
        exp(lambda_2 * t)
      C_N_t[t + 1] <- C_N0 + G * ((beta * Q_prime + lambda_2 * V_N * (beta + Q_prime)) / (beta * Q_prime * V_N * (lambda_1 - lambda_2))) * 
        exp(lambda_1 * t) - G * 
        ((beta * Q_prime + lambda_1 * V_N * (beta + Q_prime)) / (beta * Q_prime * V_N * (lambda_1 - lambda_2))) * 
        exp(lambda_2 * t)
    } else {
      # Concentration decay equations with Q replaced by Q_prime
      C_F_t[t + 1] <- C_F_t[t_g + 1] * exp(lambda_1 * (t - t_g))
      C_N_t[t + 1] <- C_N_t[t_g + 1] * exp(lambda_1 * (t - t_g))
    }
  }
  
  return(list(Time = time_points, C_F_t = C_F_t, C_N_t = C_N_t))
}

# Run Model 203
results_203 <- model_203(G, Q, Q_R, beta, gamma, epsilon_RF, epsilon_L, t_g, T, V, V_N)

combined_df <- data.frame(
    Time = rep(time_points, 4),
    Concentration = c(results_201$C_F_t, results_201$C_N_t, results_203$C_F_t, results_203$C_N_t), # 假设results_203包含了C_F_t和C_N_t
    Model = factor(rep(c("Model 201 Far Field", "Model 201 Near Field",
                         "Model 203 Far Field", "Model 203 Near Field"), each = T + 1))
)

#Comparison

ggplot(data = combined_df, aes(x = Time, y = Concentration, color = Model)) +
    geom_line() +
    labs(title = 'Comparison of Concentration vs. Time for Models 201 and 203', 
         x = 'Time (minutes)', y = 'Concentration (mg/m^3)') +
    theme_minimal() +
    scale_color_manual(values = c("Model 201 Far Field" = "blue",
                                  "Model 201 Near Field" = "red",
                                  "Model 203 Far Field" = "green",
                                  "Model 203 Near Field" = "orange")) +
    guides(color = guide_legend(title = "Model and Field"))


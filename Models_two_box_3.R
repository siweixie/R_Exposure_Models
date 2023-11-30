# Model 208

model_208 <- function(G, Q, Q_L, beta, gamma, epsilon_L, epsilon_LF) {
  epsilon_N <- Q_L / (Q_L + beta)
  
  C_F <- (gamma * G * ((1 - epsilon_L * epsilon_LF - epsilon_N * epsilon_LF * (1 - epsilon_L)))) / (Q + epsilon_LF * Q_L)
  C_N <- C_F + ((gamma * G * (1 - epsilon_L) * epsilon_N) / Q_L)
  C_L_E <- C_N + ((gamma * G * epsilon_L) / Q_L)
  C_L_F <- C_L_E * (1-epsilon_LF)
  
  return(list(C_F = C_F, C_N = C_N, C_L_E = C_L_E, C_L_F = C_L_F))
}

# Run Model 208
results_208 <- model_208(G, Q, Q_L, beta, gamma, epsilon_L, epsilon_LF)
                    

# Model 209

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

beta_i = Q_L + beta
a = V_F * V_N
b = (beta_i + Q) * V_N + beta_i * V_F
c = beta_i * (beta_i + Q) - beta_i * (beta + Q_L * (1 - epsilon_LF))
r1 = (-b + sqrt(b^2 - 4 * a*c)) / (2 * a)
r2 = (-b - sqrt(b^2 - 4 * a*c)) / (2 * a)


# Initial concentrations for the decay phase
epsilon_N <- Q_L / (Q_L + beta)
C_F0 <- (G * ((1 - epsilon_L * epsilon_LF - epsilon_N * epsilon_LF * (1 - epsilon_L)))) / (Q + epsilon_LF * Q_L)
C_N0 <- C_F0 + ((G * (1 - epsilon_L) * epsilon_N) / Q_L)


# Concentration calculations
C_F_t <- numeric(length(time_points))
C_N_t <- numeric(length(time_points))

for (t in time_points) {
  if (t <= t_g) {
    C_F_t[t + 1] <- C_F0 + ((beta_i * C_F0 - C_N0 * (beta_i + V_N * r2)) / (V_N * (r2 - r1))) *
      ((beta_i + V_N * r1) / beta_i) * exp(r1 * t) + 
      ((C_N0 * (beta_i + V_N * r1) - beta_i * C_F0) / V_N * (r2 - r1)) * 
      ((beta_i + V_N * r2) / beta_i) * exp(r2 * t)
    
    C_N_t[t + 1] <- C_N0 + ((beta_i * C_F0 - C_N0 * (beta_i + V_N * r2)) / (V_N * (r2 - r1))) *
      exp(r1 * t) +
      ((C_N0 * (beta_i + V_N * r1) - beta_i * C_F0) / V_N * (r2 - r1)) * 
      exp(r2 * t)
  } else {
    C_F_t[t + 1] <- C_F_t[t_g + 1] * exp(r1 * (t - t_g))
    C_N_t[t + 1] <- C_N_t[t_g + 1] * exp(r1 * (t - t_g))
  }
}


results_209 <- data.frame(Time = time_points, C_F_t = C_F_t, C_N_t = C_N_t)                  

                    
# Model 210

model_210 <- function(G, Q, Q_L, beta, gamma, epsilon_L, epsilon_LF, epsilon_RF, Q_R) {
  epsilon_N <- Q_L / (Q_L + beta)
  
  C_F <- (gamma * G * ((1 - epsilon_L * epsilon_LF - epsilon_L * epsilon_LF * (1 - epsilon_N)))) / 
  (Q + epsilon_RF * Q_R + epsilon_LF * Q_L)
  C_N <- C_F + ((gamma * G * (1 - epsilon_L) * epsilon_N) / Q_L)
  
  return(list(C_F = C_F, C_N = C_N))
}

# Run Model 210
results_210 <- model_210(G, Q, Q_L, beta, gamma, epsilon_L, epsilon_LF, epsilon_RF, Q_R)


# Model 211

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
Q_prime <- Q + epsilon_RF * Q_R
beta_i = Q_L + beta
r1 = (-b + sqrt(b^2 - 4 * a*c)) / (2 * a)
r2 = (-b - sqrt(b^2 - 4 * a*c)) / (2 * a)
a = V_F * V_N
b = (beta_i + Q) * V_N + beta_i * V_F
c = beta_i * (beta_i + Q_prime) - beta_i * (beta + Q_L * (1 - epsilon_LF))

# Initial concentrations for the decay phase
epsilon_N <- Q_L / (Q_L + beta)
C_F0 <- (G * ((1 - epsilon_L * epsilon_LF - epsilon_N * epsilon_LF * (1 - epsilon_L)))) / (Q_prime + epsilon_LF * Q_L)
C_N0 <- C_F0 + ((G * (1 - epsilon_L) * epsilon_N) / Q_L)


# Concentration calculations
C_F_t <- numeric(length(time_points))
C_N_t <- numeric(length(time_points))

for (t in time_points) {
  if (t <= t_g) {
    C_F_t[t + 1] <- C_F0 + ((beta_i * C_F0 - C_N0 * (beta_i + V_N * r2)) / (V_N * (r2 - r1))) *
      ((beta_i + V_N * r1) / beta_i) * exp(r1 * t) + 
      ((C_N0 * (beta_i + V_N * r1) - beta_i * C_F0) / V_N * (r2 - r1)) * 
      ((beta_i + V_N * r2) / beta_i) * exp(r2 * t)
    
    C_N_t[t + 1] <- C_N0 + ((beta_i * C_F0 - C_N0 * (beta_i + V_N * r2)) / (V_N * (r2 - r1))) *
      exp(r1 * t) +
      ((C_N0 * (beta_i + V_N * r1) - beta_i * C_F0) / V_N * (r2 - r1)) * 
      exp(r2 * t)
  } else {
    C_F_t[t + 1] <- C_F_t[t_g + 1] * exp(r1 * (t - t_g))
    C_N_t[t + 1] <- C_N_t[t_g + 1] * exp(r1 * (t - t_g))
  }
}


results_211 <- data.frame(Time = time_points, C_F_t = C_F_t, C_N_t = C_N_t)                 
                    


combined_df <- data.frame(
    Time = rep(time_points, 4),
    Concentration = c(results_209$C_F_t, results_209$C_N_t, results_211$C_F_t, results_211$C_N_t), 
    Model = factor(rep(c("Model 209 Far Field", "Model 209 Near Field",
                         "Model 211 Far Field", "Model 211 Near Field"), each = T + 1))
)

#Comparison

ggplot(data = combined_df, aes(x = Time, y = Concentration, color = Model)) +
    geom_line() +
    labs(title = 'Comparison of Concentration for Models 209 and 211', 
         x = 'Time (minutes)', y = 'Concentration (mg/m^3)') +
    theme_minimal() +
    scale_color_manual(values = c("Model 209 Far Field" = "blue",
                                  "Model 209 Near Field" = "red",
                                  "Model 211 Far Field" = "green",
                                  "Model 211 Near Field" = "orange")) +
    guides(color = guide_legend(title = "Model and Field"))          


near_field_df <- data.frame(
    Time = rep(time_points, 6),  
    Concentration = c(results_201$C_N_t, results_203$C_N_t, results_205$C_N_t, 
                      results_207$C_N_t, results_209$C_N_t, results_211$C_N_t), 
    Model = factor(rep(c("Model 201 Near Field", "Model 203 Near Field",
                         "Model 205 Near Field", "Model 207 Near Field",
                         "Model 209 Near Field", "Model 211 Near Field"), each = T + 1))
)

far_field_df <- data.frame(
    Time = rep(time_points, 6),  
    Concentration = c(results_201$C_F_t, results_203$C_F_t, results_205$C_F_t, 
                      results_207$C_F_t, results_209$C_F_t, results_211$C_F_t), 
    Model = factor(rep(c("Model 201 Far Field", "Model 203 Far Field",
                         "Model 205 Far Field", "Model 207 Far Field",
                         "Model 209 Far Field", "Model 211 Far Field"), each = T + 1))

  
# Plot for Near Field Data
ggplot(data = near_field_df, aes(x = Time, y = Concentration, color = Model)) +
    geom_line() +
    labs(title = 'Concentration Comparison for Models - Near Field', 
         x = 'Time (minutes)', y = 'Concentration (mg/m^3)') +
    theme_minimal() +
    scale_color_brewer(palette = "Set1") +
    guides(color = guide_legend(title = "Model"))

# Plot for Far Field Data
ggplot(data = far_field_df, aes(x = Time, y = Concentration, color = Model)) +
    geom_line() +
    labs(title = 'Concentration Comparison for Models - Far Field', 
         x = 'Time (minutes)', y = 'Concentration (mg/m^3)') +
    theme_minimal() +
    scale_color_brewer(palette = "Set2") +
    guides(color = guide_legend(title = "Model"))

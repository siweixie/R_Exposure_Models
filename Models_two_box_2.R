# Model 204

model_204 <- function(G, Q, Q_L, beta, gamma, epsilon_L) {
  epsilon_N <- Q_L / (Q_L + beta)
  
  C_F <- (gamma * G * (1 - epsilon_L) * (1 - epsilon_N)) / (Q + Q_L)
  C_N <- C_F + ((gamma * G * (1 - epsilon_L) * epsilon_N) / Q_L)
  C_L_E <- C_N + ((gamma * G * epsilon_L) / Q_L)
  
  return(list(C_F = C_F, C_N = C_N, C_L_E = C_L_E))
}

# Run Model 204
results_204 <- model_204(G, Q, Q_L, beta, gamma, epsilon_L)


# Define model 205

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
c = beta_i * (beta_i + Q) - beta * beta_i
r1 = (-b + sqrt(b^2 - 4 * a*c)) / (2 * a)
r2 = (-b - sqrt(b^2 - 4 * a*c)) / (2 * a)

# Initial concentrations for the decay phase
epsilon_N <- Q_L / (Q_L + beta)
C_F0 <- (G * (1 - epsilon_L) * (1 - epsilon_N)) / (Q + Q_L) # Steady state far field concentration
C_N0 <- C_F0 + ((G * (1 - epsilon_L) * epsilon_N) / Q_L) # Steady state near field concentration


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


results_205 <- data.frame(Time = time_points, C_F_t = C_F_t, C_N_t = C_N_t)


# Model 206

model_206 <- function(G, Q, Q_L, epsilon_RF, Q_R, beta, gamma, epsilon_L) {
  epsilon_N <- Q_L / (Q_L + beta)
  
  C_F <- (gamma * G * (1 - epsilon_L) * (1 - epsilon_N)) / (Q + epsilon_RF * Q_R + Q_L)
  C_N <- C_F + ((gamma * G * (1 - epsilon_L) * epsilon_N) / Q_L)
  C_L_E <- C_F + ((gamma * G * epsilon_L) / Q_L)
  C_R_F <- C_F * (1 - epsilon_RF)
  
  return(list(C_F = C_F, C_N = C_N, C_L_E = C_L_E, C_R_F = C_R_F))
}

# Run Model 206
results_206 <- model_206(G, Q, Q_L, epsilon_RF, Q_R, beta, gamma, epsilon_L)


# Define model 207

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
b = (beta_i + Q_prime) * V_N + beta_i * V_F
c = beta_i * (beta_i + Q_prime) - beta * beta_i

# Initial concentrations for the decay phase
epsilon_N <- Q_L / (Q_L + beta)
C_F0 <- (G * (1 - epsilon_L) * (1 - epsilon_N)) / (Q_prime + Q_L) # Steady state far field concentration
C_N0 <- C_F0 + ((G * (1 - epsilon_L) * epsilon_N) / Q_L) # Steady state near field concentration


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


results_207 <- data.frame(Time = time_points, C_F_t = C_F_t, C_N_t = C_N_t)



combined_df <- data.frame(
    Time = rep(time_points, 4),
    Concentration = c(results_205$C_F_t, results_205$C_N_t, results_207$C_F_t, results_207$C_N_t), 
    Model = factor(rep(c("Model 205 Far Field", "Model 205 Near Field",
                         "Model 207 Far Field", "Model 207 Near Field"), each = T + 1))
)

#Comparison

ggplot(data = combined_df, aes(x = Time, y = Concentration, color = Model)) +
    geom_line() +
    labs(title = 'Comparison of Concentration for Models 205 and 207', 
         x = 'Time (minutes)', y = 'Concentration (mg/m^3)') +
    theme_minimal() +
    scale_color_manual(values = c("Model 205 Far Field" = "blue",
                                  "Model 205 Near Field" = "red",
                                  "Model 207 Far Field" = "green",
                                  "Model 207 Near Field" = "orange")) +
    guides(color = guide_legend(title = "Model and Field"))


# MSE and R-squared for Near Field
mse_205_207_near <- mean((results_205$C_N_t - results_207$C_N_t)^2)
r2_205_207_near <- cor(results_205$C_N_t, results_207$C_N_t)^2

# Normality Test for Near Field
shapiro_205_near <- shapiro.test(results_205$C_N_t)
shapiro_207_near <- shapiro.test(results_207$C_N_t)

# Wilcoxon Test for Near Field (if needed)
if (shapiro_205_near$p.value < 0.05 || shapiro_207_near$p.value < 0.05) {
    wilcox_205_207_near <- wilcox.test(results_205$C_N_t, results_207$C_N_t, paired = TRUE)
}

# Output results
mse_205_207_near
r2_205_207_near
shapiro_205_near
shapiro_207_near
if (exists("wilcox_205_207_near")) wilcox_205_207_near


# MSE and R-squared for Far Field
mse_205_207_far <- mean((results_205$C_F_t - results_207$C_F_t)^2)
r2_205_207_far <- cor(results_205$C_F_t, results_207$C_F_t)^2

# Normality Test for Far Field
shapiro_205_far <- shapiro.test(results_205$C_F_t)
shapiro_207_far <- shapiro.test(results_207$C_F_t)

# Wilcoxon Test for Far Field (if needed)
if (shapiro_205_far$p.value < 0.05 || shapiro_207_far$p.value < 0.05) {
    wilcox_205_207_far <- wilcox.test(results_205$C_F_t, results_207$C_F_t, paired = TRUE)
}

# Output results
mse_205_207_far
r2_205_207_far
shapiro_205_far
shapiro_207_far
if (exists("wilcox_205_207_far")) wilcox_205_207_far


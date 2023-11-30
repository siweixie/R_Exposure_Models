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
QL <- 5     # m^3/min
epsilon_L <- 0.5  # Local control efficiency
epsilon_LF <- 0.75 # Filtration efficiency
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

# Data preparation for plotting
df <- data.frame(Time = time_points, Far_Field = C_F_t, Near_Field = C_N_t)

# Plot the results
ggplot(data = df, aes(x = Time)) +
  geom_line(aes(y = Far_Field, color = "Far Field")) +
  geom_line(aes(y = Near_Field, color = "Near Field")) +
  labs(title = 'Concentration vs. Time for a Cyclic Task', 
       x = 'Time (minutes)', y = 'Concentration (mg/m^3)') +
  theme_minimal() +
  scale_color_manual("", 
                     breaks = c("Far Field", "Near Field"),
                     values = c("blue", "red"))


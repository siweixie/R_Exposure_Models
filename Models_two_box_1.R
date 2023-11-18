# Define model 200
model_200 <- function(G, Q, beta, gamma) {
  C_F <- gamma * G / Q # Far field
  C_N <- C_F + ((gamma * G) / beta) # Near field
  return(list(C_F = C_F, C_N = C_N))
}

# Define model 201
model_201 <- function(G, Q, V, V_N, beta, t_g, T) {
  # Create a time vector
  time_vector <- seq(0, T, by = 1)

  # Lambda calculations
  lambda_1 <- 0.5 * (-(beta * V_F + V_N * (beta + Q)) / (V_N * V_F)) + sqrt(((beta * V_F + V_N * (beta + Q)) / (V_N * V_F))^2 - 4 * ((beta * Q) / (V_N * V_F)))
  lambda_2 <- 0.5 * (-(beta * V_F + V_N * (beta + Q)) / (V_N * V_F)) - sqrt(((beta * V_F + V_N * (beta + Q)) / (V_N * V_F))^2 - 4 * ((beta * Q) / (V_N * V_F)))
  
  C_F_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return(
        C_F + G * ((lambda_1 * V_N + beta) / beta) * 
        ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
        exp(lambda_1 * t) - G * ((lambda_2 * V_N + beta) / beta) *
        ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) *
        exp(lambda_2 * t)
      )
    } else {
      return(NA)
    }
  })
  
  C_F0 <- (G * (1 - epsilon_L)) / (Q + Q_L) * (1 - exp((-(Q + Q_L) * t_g) / V))
  
  C_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(C0 * exp((-(Q + Q_L) * (t - t_g)) / V))
    } else {
      return(NA)
    }
  })
  
  return(list(time = time_vector, concentration_rise = C_rise, concentration_decay = C_decay))
}




# Model 201 equations for concentration rise (G > 0)
C_F_t <- function(t) {
  C_F + G * ((lambda_1 * V_N + beta) / beta) * 
  ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
  exp(lambda_1 * t) - G * ((lambda_2 * V_N + beta) / beta) *
  ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) *
  exp(lambda_2 * t)
}

C_N_t <- function(t) {
  C_N + G * (beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2)) * exp(lambda_1 * t) - G * (lambda_2 * V_N + beta) / beta * exp(lambda_2 * t)
}

# Model 201 equations for concentration decay (G = 0)
C_F_t_decay <- function(t) {
  (lambda_1 * V_N + beta) / (beta * V_N * (lambda_1 - lambda_2)) * (beta * (C_F0 - C_N0) - lambda_2 * V_N * C_N0) * exp(lambda_1 * t) + 
  (lambda_2 * V_N + beta) / (beta * V_N * (lambda_1 - lambda_2)) * (beta * (C_N0 - C_F0) - lambda_1 * V_N * C_N0) * exp(lambda_2 * t)
}

C_N_t_decay <- function(t) {
  (beta * (C_F0 - C_N0) - lambda_2 * V_N * C_N0) / (V_N * (lambda_1 - lambda_2)) * exp(lambda_1 * t) + 
  (beta * (C_N0 - C_F0) - lambda_1 * V_N * C_N0) / (V_N * (lambda_1 - lambda_2)) * exp(lambda_2 * t)
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

results_200 <- model_200(G, Q, beta, gamma)


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

   C_N_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return(
        C_F + G * ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
        exp(lambda_1 * t) - G * 
        ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) *
        exp(lambda_2 * t)
      )
    } else {
      return(NA)
    }
  })
  
  
  C_F0 <-  C_F + G * ((lambda_1 * V_N + beta) / beta) * 
           ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
           exp(lambda_1 * t_g) - G * ((lambda_2 * V_N + beta) / beta) *
           ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) *
           exp(lambda_2 * t_g)

  C_N0 <- C_F + G * ((beta * Q + lambda_2 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) * 
          exp(lambda_1 * t) - G * 
          ((beta * Q + lambda_1 * V_N * (beta + Q)) / (beta * Q * V_N * (lambda_1 - lambda_2))) *
          exp(lambda_2 * t)
  
  C_F_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(
        (((lambda_1 * V_N + beta) * (beta * (C_F0 - C_N0) - (lambda_2 * V_N * C_N0)))
        / (beta * V_N * (lambda_1 - lambda_2))) * exp(lambda_1 * (t - t_g)) +
        (((lambda_2 * V_N + beta) * (beta * (C_F0 - C_N0) - (lambda_1 * V_N * C_N0)))
        / (beta * V_N * (lambda_1 - lambda_2))) * exp(lambda_2 * (t - t_g))
      )
    } else {
      return(NA)
    }
  })

  C_N_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(
        (((beta * (C_F0 - C_N0) - (lambda_2 * V_N * C_N0)))
        / (V_N * (lambda_1 - lambda_2))) * exp(lambda_1 * (t - t_g)) +
        (((beta * (C_F0 - C_N0) - (lambda_1 * V_N * C_N0)))
        / (V_N * (lambda_1 - lambda_2))) * exp(lambda_2 * (t - t_g))
      )
    } else {
      return(NA)
    }
  })
  
  return(list(time = time_vector, concentration_F_rise = C_F_rise, concentration_N_rise = C_N_rise, 
              concentration_F_decay = C_F_decay, concentration_N_decay = C_N_decay))
}

# Define the model parameters with units
t <- 60    # minutes
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

results_200 <- model_200(G, Q, beta, gamma)



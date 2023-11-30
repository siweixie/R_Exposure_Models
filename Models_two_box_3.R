# Model 208

model_208 <- function(G, Q, Q_L, beta, gamma, epsilon_L, epsilon_LF) {
  epsilon_N <- Q_L / (Q_L + beta)
  
  C_F <- (gamma * G * ((1 - epsilon_L * epsilon_LF - epsilon_N * epsilon_LF * (1 - epsilon_L)))) / (Q + epsilon_LF * Q_L)
  C_N <- C_F + ((gamma * G * (1 - epsilon_L) * epsilon_N) / Q_L)
  C_L_E <- C_N + ((gamma * G * epsilon_L) / Q_L)
  C_L_F <- C_L_E * (1-epsilon_LF
  
  return(list(C_F = C_F, C_N = C_N, C_L_E = C_L_E, C_L_F = C_L_F))
}

# Run Model 208
results_208 <- model_208(G, Q, Q_L, beta, gamma, epsilon_L, epsilon_LF)
                    

# Model 209
                   
                    
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


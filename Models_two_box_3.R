# Model 208

model_208 <- function(G, Q, Q_L, beta, gamma, epsilon_L) {
  epsilon_N <- Q_L / (Q_L + beta)
  
  C_F <- (gamma * G * (1 - epsilon_L) * (1 - epsilon_N)) / (Q + Q_L)
  C_N <- C_F + ((gamma * G * (1 - epsilon_L) * epsilon_N) / Q_L)
  C_L_E <- C_N + ((gamma * G * epsilon_L) / Q_L)
  
  return(list(C_F = C_F, C_N = C_N, C_L_E = C_L_E))
}

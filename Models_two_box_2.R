# Model 204

model_204 <- function(G, Q, QL, beta, gamma, epsilon_L) {
  epsilon_N <- QL / (QL + beta)
  
  C_F <- gamma * G * (1 - epsilon_L) * (1 - epsilon_N) / (Q + QL)
  C_N <- C_F + gamma * G * epsilon_L * epsilon_N / QL
  C_L_E <- C_N + gamma * G * epsilon_L / QL
  
  return(list(C_F = C_F, C_N = C_N, C_L_E = C_L_E))
}

# Run Model 204
results_204 <- model_204(G, Q, QL, beta, gamma, epsilon_L)

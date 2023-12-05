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
  return(model_101(G, Q_instead, V, t_g, T))
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

# Monte Carlo

G_range <- c(60, 150)
Q_range <- c(10, 30)
Q_R_range <- c(2, 10)
V_range <- c(60, 150)
t_g_range <- c(5, 30)
T <- 60
gamma_range <- t_g_range / T
epsilon_RF_range <- c(0.5, 0.95)

# simulation times
n_simulations <- 1000

results_101 <- vector("list", n_simulations)
results_103 <- vector("list", n_simulations)

# Monte Carlo simulation
set.seed(123)
for (i in 1:n_simulations) {
    # Randomly selection
    G <- runif(1, G_range[1], G_range[2])
    Q <- runif(1, Q_range[1], Q_range[2])
    Q_R <- runif(1, Q_R_range[1], Q_R_range[2])
    V <- runif(1, V_range[1], V_range[2])
    t_g <- runif(1, t_g_range[1], t_g_range[2])
    gamma <- runif(1, gamma_range[1], gamma_range[2])
    epsilon_RF <- runif(1, epsilon_RF_range[1], epsilon_RF_range[2])
    
    results_101[[i]] <- model_101(G, Q, V, t_g, T)
    results_103[[i]] <- model_103(G, Q, Q_R, epsilon_RF, V, t_g, T)
}

# Transfer to data.frame
df_101 <- do.call(rbind, lapply(results_101, function(x) data.frame(Time = x$time, Concentration = c(x$concentration_rise, x$concentration_decay))))
df_103 <- do.call(rbind, lapply(results_103, function(x) data.frame(Time = x$time, Concentration = c(x$concentration_rise, x$concentration_decay))))

                                
average_and_ci <- function(df, n) {
    df_summary <- aggregate(Concentration ~ Time, data = df, mean)
    df_summary$SD <- aggregate(Concentration ~ Time, data = df, sd)$Concentration
    df_summary$CI_lower <- df_summary$Concentration - qt(0.975, df=n-1) * df_summary$SD / sqrt(n)
    df_summary$CI_upper <- df_summary$Concentration + qt(0.975, df=n-1) * df_summary$SD / sqrt(n)
    return(df_summary)
}

average_101 <- average_and_ci(df_101, n_simulations)
average_103 <- average_and_ci(df_103, n_simulations)

# 绘制平均曲线和置信区间
ggplot() +
    geom_ribbon(data = average_101, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "red", alpha = 0.2) +
    geom_line(data = average_101, aes(x = Time, y = Concentration), color = "red") +
    geom_ribbon(data = average_103, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "blue", alpha = 0.2) +
    geom_line(data = average_103, aes(x = Time, y = Concentration), color = "blue") +
    labs(title = "Model 101 vs 103 with Monte Carlo Simulation",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal()






# Define Model 104
model_104 <- function(G, Q, Q_L, epsilon_L, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L)) / (Q + Q_L)
  C_L_avg <- C_avg + (epsilon_L * gamma * G) / Q_L
  return(list(C_avg = C_avg, C_L_avg = C_L_avg))
}


# Define Model 105
model_105 <- function(G, Q, Q_L, epsilon_L, V, t_g, T) {
  # Create a time vector
  time_vector <- seq(0, T, by = 1)
  
  C_rise <- sapply(time_vector, function(t) {
    if (t <= t_g) {
      return(((G * (1 - epsilon_L)) / (Q + Q_L)) * (1 - exp((-(Q + Q_L) * t) / V)))
    } else {
      return(NA)
    }
  })
  
  C0 <- (G * (1 - epsilon_L)) / (Q + Q_L) * (1 - exp((-(Q + Q_L) * t_g) / V))
  
  C_decay <- sapply(time_vector, function(t) {
    if (t > t_g) {
      return(C0 * exp((-(Q + Q_L) * (t - t_g)) / V))
    } else {
      return(NA)
    }
  })
  
  return(list(time = time_vector, concentration_rise = C_rise, concentration_decay = C_decay))
}


# Define Model 106
model_106 <- function(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, gamma) {
  C_avg <- (gamma * G * (1 - epsilon_L)) / ((Q + epsilon_RF * Q_R) + Q_L)
  C_L_avg <- C_avg + (epsilon_L * gamma * G) / Q_L
  C_RF <- C_avg * (1 - epsilon_RF)
  return(list(C_avg = C_avg, C_L_avg = C_L_avg, C_RF = C_RF))
}


# Define Model 107
model_107 <- function(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, V, t_g, T) {
  Q_instead <- Q + epsilon_RF * Q_R
  return(model_105(G, Q_instead, Q_L, epsilon_L, V, t_g, T))
}


T <- 60   # Total time (minutes)
t_g <- 15 # Time of generation (minutes)
G <- 100  # mg/min
Q <- 20   # m^3/min
Q_L <- 5  # m^3/min
epsilon_L <- 0.5  # Efficiency of local exhaust
epsilon_L_F <- 0.75  # Efficiency of local exhaust filtration
Q_R <- 5  # m^3/min
epsilon_RF <- 0.9  # Efficiency of recirculation filtration
V <- 100  # m^3
gamma <- 0.25 


# Monte Carlo

G_range <- c(60, 150)
Q_range <- c(10, 30)
Q_L_range <- c(2, 8)
epsilon_L_range <- c(0.2, 0.8)
epsilon_L_F_range <- c(0.5, 0.95)
Q_R_range <- c(2, 10)
epsilon_RF_range <- c(0.5, 0.95)
V_range <- c(60, 150)
t_g_range <- c(5, 30)
T <- 60
gamma_range <- t_g_range / T


# simulation times
n_simulations <- 1000

results_105 <- vector("list", n_simulations)
results_107 <- vector("list", n_simulations)

# Monte Carlo simulation
set.seed(123)
for (i in 1:n_simulations) {
    # Randomly selection
    G <- runif(1, G_range[1], G_range[2])
    Q <- runif(1, Q_range[1], Q_range[2])
    Q_L <- runif(1, Q_L_range[1], Q_L_range[2])
    epsilon_L <- runif(1, epsilon_L_range[1], epsilon_L_range[2])
    epsilon_L_F <- runif(1, epsilon_L_F_range[1], epsilon_L_F_range[2])
    Q_R <- runif(1, Q_R_range[1], Q_R_range[2])
    V <- runif(1, V_range[1], V_range[2])
    t_g <- runif(1, t_g_range[1], t_g_range[2])
    gamma <- runif(1, gamma_range[1], gamma_range[2])
    epsilon_RF <- runif(1, epsilon_RF_range[1], epsilon_RF_range[2])
    
    results_105[[i]] <- model_105(G, Q, Q_L, epsilon_L, V, t_g, T)
    results_107[[i]] <- model_107(G, Q, Q_L, epsilon_L, Q_R, epsilon_RF, V, t_g, T)
}

# Transfer to data.frame
df_105 <- do.call(rbind, lapply(results_105, function(x) data.frame(Time = x$time, 
Concentration = c(x$concentration_rise, x$concentration_decay))))

df_107 <- do.call(rbind, lapply(results_107, function(x) data.frame(Time = x$time, 
Concentration = c(x$concentration_rise, x$concentration_decay))))

                                
average_and_ci <- function(df, n) {
    df_summary <- aggregate(Concentration ~ Time, data = df, mean)
    df_summary$SD <- aggregate(Concentration ~ Time, data = df, sd)$Concentration
    df_summary$CI_lower <- df_summary$Concentration - qt(0.975, df=n-1) * df_summary$SD / sqrt(n)
    df_summary$CI_upper <- df_summary$Concentration + qt(0.975, df=n-1) * df_summary$SD / sqrt(n)
    return(df_summary)
}

average_105 <- average_and_ci(df_105, n_simulations)
average_107 <- average_and_ci(df_107, n_simulations)


ggplot() +
    geom_ribbon(data = average_105, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "red", alpha = 0.2) +
    geom_line(data = average_105, aes(x = Time, y = Concentration), color = "red") +
    geom_ribbon(data = average_107, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "blue", alpha = 0.2) +
    geom_line(data = average_107, aes(x = Time, y = Concentration), color = "blue") +
    labs(title = "Model 105 vs 107 with Monte Carlo Simulation",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal()
  theme_minimal()


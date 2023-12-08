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

#Visualization
library(ggplot2)

# Monte Carlo

# simulation times
n_simulations <- 1000

results_101 <- vector("list", n_simulations)
results_103 <- vector("list", n_simulations)

# Monte Carlo simulation
set.seed(123)  

for (i in 1:n_simulations) {

    G <- rlnorm(1, log(26000), sqrt(log(2)))
    

    Q <- rnorm(1, 205, 6.2)
    Q_L <- rnorm(1, 25, 0.75)
    

    eL_min <- 0.2 - (1.96 * 0.4)
    eL_max <- 0.2 + (1.96 * 0.4)
    eLF_min <- 0.5 - (1.96 * 0.9)
    eLF_max <- 0.5 + (1.96 * 0.9)
    eRF_min <- 0.5 - (1.96 * 0.9)
    eRF_max <- 0.5 + (1.96 * 0.9)
    

    epsilon_L <- runif(1, eL_min, eL_max)
    epsilon_L_F <- runif(1, eLF_min, eLF_max)
    epsilon_RF <- runif(1, eRF_min, eRF_max)
    

    V <- 975  
    t_g <- 120  
    T <- 480  
    Q_R <- 20 + rnorm(1, 0, 0.6)  
    
    results_101[[i]] <- model_101(G, Q, V, t_g, T)
    results_103[[i]] <- model_103(G, Q, Q_R, epsilon_RF, V, t_g, T)
}

# Transfer to data.frame
df_101 <- do.call(rbind, lapply(results_101, function(x) data.frame(Time = x$time, 
                                                                    Concentration = c(x$concentration_rise, x$concentration_decay))))
df_103 <- do.call(rbind, lapply(results_103, function(x) data.frame(Time = x$time, 
                                                                    Concentration = c(x$concentration_rise, x$concentration_decay))))


average_and_ci <- function(df, n) {
    df_summary <- aggregate(Concentration ~ Time, data = df, mean)
    df_summary$SD <- aggregate(Concentration ~ Time, data = df, sd)$Concentration
    df_summary$CI_lower <- df_summary$Concentration - qt(0.975, df=n-1) * df_summary$SD / sqrt(n)
    df_summary$CI_upper <- df_summary$Concentration + qt(0.975, df=n-1) * df_summary$SD / sqrt(n)
    return(df_summary)
}

average_101 <- average_and_ci(df_101, n_simulations)
average_103 <- average_and_ci(df_103, n_simulations)

ggplot() +
    geom_ribbon(data = average_101, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "red", alpha = 0.2) +
    geom_line(data = average_101, aes(x = Time, y = Concentration), color = "red") +
    geom_ribbon(data = average_103, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "blue", alpha = 0.2) +
    geom_line(data = average_103, aes(x = Time, y = Concentration), color = "blue") +
    labs(title = "Model 101 vs 103 with Monte Carlo Simulation",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal()  

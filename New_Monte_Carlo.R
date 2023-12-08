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
    Q_R <- rnorm(1, 200, 6)  
    
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


results_105 <- vector("list", n_simulations)
results_107 <- vector("list", n_simulations)

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
    Q_R <- rnorm(1, 200, 6)  
    
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



# Define Model 108
model_108 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, gamma) {
    C_avg <- (gamma * G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)
    C_LE <- C_avg + (epsilon_L * gamma * G) / Q_L
    C_LF <- C_LE * (1 - epsilon_L_F)
    return(list(C_avg = C_avg, C_LE = C_LE, C_LF = C_LF))
}


# Define Model 109
model_109 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, V, t_g, T) {
    # Create a time vector
    time_vector <- seq(0, T, by = 1)
    
    C_rise <- sapply(time_vector, function(t) {
        if (t <= t_g) {
            return(((G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)) * (1 - exp((-(Q + epsilon_L_F * Q_L) * t) / V)))
        } else {
            return(NA)
        }
    })
    
    C0 <- ((G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)) * (1 - exp((-(Q + epsilon_L_F * Q_L) * t_g) / V))
    
    C_decay <- sapply(time_vector, function(t) {
        if (t > t_g) {
            return(C0 * exp((-(Q + epsilon_L_F * Q_L) * (t - t_g)) / V))
        } else {
            return(NA)
        }
    })
    
    return(list(time = time_vector, concentration_rise = C_rise, concentration_decay = C_decay))
}


# Define Model 110
model_110 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, gamma) {
    C_avg <- (gamma * G * (1 - epsilon_L * epsilon_L_F)) / ((Q + epsilon_RF * Q_R) + epsilon_L_F * Q_L)
    C_LE <- C_avg + (epsilon_L * gamma * G) / Q_L
    C_LF <- C_LE * (1 - epsilon_L_F)
    C_RF <- C_avg * (1 - epsilon_RF)
    return(list(C_avg = C_avg, C_LE = C_LE, C_LF = C_LF, C_RF = C_RF))
}


# Define Model 111
model_111 <- function(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, V, t_g, T) {
    Q_instead <- Q + epsilon_RF * Q_R
    return(model_109(G, Q_instead, Q_L, epsilon_L, epsilon_L_F, V, t_g, T))
}


results_109 <- vector("list", n_simulations)
results_111 <- vector("list", n_simulations)

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
    Q_R <- rnorm(1, 200, 6)  
    
    results_109[[i]] <- model_109(G, Q, Q_L, epsilon_L, epsilon_L_F, V, t_g, T)
    results_111[[i]] <- model_111(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, V, t_g, T)
}



# Transfer to data.frame
df_109 <- do.call(rbind, lapply(results_109, function(x) data.frame(Time = x$time, 
                                                                    Concentration = c(x$concentration_rise, x$concentration_decay))))

df_111 <- do.call(rbind, lapply(results_111, function(x) data.frame(Time = x$time, 
                                                                    Concentration = c(x$concentration_rise, x$concentration_decay))))

# Calculate Average and Confidence Interval
average_109 <- average_and_ci(df_109, n_simulations)
average_111 <- average_and_ci(df_111, n_simulations)


# Plotting with ggplot
ggplot() +
    geom_ribbon(data = average_109, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "red", alpha = 0.2) +
    geom_line(data = average_109, aes(x = Time, y = Concentration), color = "red") +
    geom_ribbon(data = average_111, aes(x = Time, ymin = CI_lower, ymax = CI_upper), fill = "blue", alpha = 0.2) +
    geom_line(data = average_111, aes(x = Time, y = Concentration), color = "blue") +
    labs(title = "Model 109 vs 111 with Monte Carlo Simulation",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal()


# Overall comparison
average_101$Model <- "Model 101"
average_103$Model <- "Model 103"
average_105$Model <- "Model 105"
average_107$Model <- "Model 107"
average_109$Model <- "Model 109"
average_111$Model <- "Model 111"

combined_averages <- rbind(average_101, average_103, average_105, average_107, average_109, average_111)

ggplot(combined_averages, aes(x = Time, y = Concentration, color = Model)) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = Model), alpha = 0.2) +
    geom_line() +
    labs(title = "Overall Model Comparison with Monte Carlo Simulation",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal()  

# Or

ggplot(combined_averages, aes(x = Time, y = Concentration, group = Model)) +
    geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper, fill = Model), alpha = 0.2) +
    geom_line(aes(color = Model)) +
    scale_fill_manual(values = c("Model 101" = "red", "Model 103" = "blue", "Model 105" = "yellow", 
                                 "Model 107" = "green", "Model 109" = "purple", "Model 111" = "pink")) +
    scale_color_manual(values = c("Model 101" = "red", "Model 103" = "blue", "Model 105" = "yellow", 
                                  "Model 107" = "green", "Model 109" = "purple", "Model 111" = "pink")) +
    labs(title = "Overall Model Comparison with Monte Carlo Simulation",
         x = "Time (minutes)", 
         y = "Concentration (mg/m^3)") +
    theme_minimal() +
    guides(fill = guide_legend(title = "Model"), color = guide_legend(title = "Model"))                                

library(ggplot2)

set.seed(123)

G_values <- rlnorm(1000, log(26000), sqrt(log(2)))

Q_values <- rnorm(1000, 205, 6.2)
Q_L_values <- rnorm(1000, 25, 0.75)

eL_min <- 0.2 - (1.96 * 0.4)
eL_max <- 0.2 + (1.96 * 0.4)
eLF_min <- 0.5 - (1.96 * 0.9)
eLF_max <- 0.5 + (1.96 * 0.9)
eRF_min <- 0.5 - (1.96 * 0.9)
eRF_max <- 0.5 + (1.96 * 0.9)

epsilon_L_values <- runif(1000, eL_min, eL_max)
epsilon_L_F_values <- runif(1000, eLF_min, eLF_max)
epsilon_RF_values <- runif(1000, eRF_min, eRF_max)

QR_values <- rnorm(1000, 20, 0.6)

p1 <- ggplot(data.frame(G = G_values), aes(x = G)) +
    geom_histogram(aes(y = ..density..), binwidth = 200, fill = "blue", alpha = 0.5) +
    geom_density(color = "red", adjust = 1) +
    labs(title = "Distribution of G", x = "G values", y = "Density")

p2 <- ggplot(data.frame(Q = Q_values), aes(x = Q)) +
    geom_histogram(aes(y = ..density..), binwidth = 2, fill = "blue", alpha = 0.5) +
    geom_density(color = "red", adjust = 1) +
    labs(title = "Distribution of Q", x = "Q values", y = "Density")

p3 <- ggplot(data.frame(Q_L = Q_L_values), aes(x = Q_L)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.1, fill = "blue", alpha = 0.5) +
    geom_density(color = "red", adjust = 1) +
    labs(title = "Distribution of Q_L", x = "Q_L values", y = "Density")

p4 <- ggplot(data.frame(epsilon_L = epsilon_L_values), aes(x = epsilon_L)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.01, fill = "blue", alpha = 0.5) +
    labs(title = "Distribution of epsilon_L", x = "epsilon_L values", y = "Density")

p5 <- ggplot(data.frame(epsilon_L_F = epsilon_L_F_values), aes(x = epsilon_L_F)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.02, fill = "green", alpha = 0.5) +
    labs(title = "Distribution of epsilon_L_F", x = "epsilon_L_F values", y = "Density")

p6 <- ggplot(data.frame(epsilon_RF = epsilon_RF_values), aes(x = epsilon_RF)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.02, fill = "purple", alpha = 0.5) +
    labs(title = "Distribution of epsilon_RF", x = "epsilon_RF values", y = "Density")

p7 <- ggplot(data.frame(QR = QR_values), aes(x = QR)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.1, fill = "orange", alpha = 0.5) +
    geom_density(color = "red", adjust = 1) +
    labs(title = "Distribution of QR", x = "QR values", y = "Density")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)

# Define model 200
C_F <- gamma * G / Q # Far field
C_N <- C_F + ((gamma * G) / beta) # Near field

# Define model 201


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

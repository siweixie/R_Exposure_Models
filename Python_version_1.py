import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import shapiro, wilcoxon

# Define Model 100 with the gamma factor
def model_100(G, Q, gamma):
    C_steady = gamma * G / Q
    return C_steady

# Define Model 101
def model_101(G, Q, V, t_g, T):
    time_vector = np.arange(0, T+1)
    C_rise = [(G / Q) * (1 - np.exp((-Q * t) / V)) if t <= t_g else np.nan for t in time_vector]
    C0 = (G / Q) * (1 - np.exp(-Q * t_g / V))
    C_decay = [C0 * np.exp(-Q * (t - t_g) / V) if t > t_g else np.nan for t in time_vector]
    return time_vector, np.array(C_rise), np.array(C_decay)

# Define Model 102
def model_102(G, Q, Q_R, epsilon_RF, gamma):
    C_bar = (gamma * G) / (Q + epsilon_RF * Q_R)
    C_RF = C_bar * (1 - epsilon_RF)
    return C_bar, C_RF

# Define Model 103
def model_103(G, Q, Q_R, epsilon_RF, V, t_g, T):
    Q_instead = Q + epsilon_RF * Q_R
    return model_101(G, Q_instead, V, t_g, T)

# Set the parameters
T = 60   # Total time (minutes)
t_g = 15 # Time of generation (minutes)
G = 100  # mg/min
Q = 20   # m^3/min
epsilon_RF = 0.9  # Efficiency of recirculation filtration
V = 100  # m^3
gamma = 0.25  
Q_R = 5  # m^3/min

# Perform calculations for models
C_steady_100 = model_100(G, Q, gamma)
time_101, C_rise_101, C_decay_101 = model_101(G, Q, V, t_g, T)
C_bar_102, C_RF_102 = model_102(G, Q, Q_R, epsilon_RF, gamma)
time_103, C_rise_103, C_decay_103 = model_103(G, Q, Q_R, epsilon_RF, V, t_g, T)

# Remove NaN values for statistical tests and visualization
valid_indices_101 = ~np.isnan(C_rise_101)
valid_indices_103 = ~np.isnan(C_rise_103)

# Visualization: Model 101 vs Model 103
plt.figure(figsize=(10, 5))
plt.plot(time_101[valid_indices_101], C_rise_101[valid_indices_101], label='Model 101 Rise')
plt.plot(time_101[~valid_indices_101], C_decay_101[~valid_indices_101], label='Model 101 Decay')
plt.plot(time_103[valid_indices_103], C_rise_103[valid_indices_103], label='Model 103 Rise', linestyle='--')
plt.plot(time_103[~valid_indices_103], C_decay_103[~valid_indices_103], label='Model 103 Decay', linestyle='--')
plt.axhline(y=C_steady_100, color='red', linestyle='-.', label='Model 100 Steady State')
plt.title('Model 101 vs Model 103')
plt.xlabel('Time (minutes)')
plt.ylabel('Concentration (mg/m^3)')
plt.legend()
plt.show()

# Statistical Analysis
# Remove NaN values for statistical tests
C_101 = C_rise_101[valid_indices_101]
C_103 = C_rise_103[valid_indices_103]

# Normality test
print('Shapiro test for Model 101:', shapiro(C_101))
print('Shapiro test for Model 103:', shapiro(C_103))

# If it is not a normal distribution (p < 0.05), using Wilcoxon test
p_value_101 = shapiro(C_101).pvalue
p_value_103 = shapiro(C_103).pvalue
if p_value_101 < 0.05 and p_value_103 < 0.05:
    w_stat, p_value_wilcox = wilcoxon(C_101, C_103)
    print('Wilcoxon test statistic:', w_stat)
    print('Wilcoxon test p-value:', p_value_wilcox)


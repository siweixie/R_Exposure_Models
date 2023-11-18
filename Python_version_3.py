import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import shapiro, wilcoxon

# Define the models in Python
def model_108(G, Q, Q_L, epsilon_L, epsilon_L_F, gamma):
    C_avg = (gamma * G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L)
    C_LE = C_avg + (epsilon_L * gamma * G) / Q_L
    C_LF = C_LE * (1 - epsilon_L_F)
    return C_avg, C_LE, C_LF

def model_109(G, Q, Q_L, epsilon_L, epsilon_L_F, V, t_g, T):
    time_vector = np.arange(0, T + 1)
    C_rise = [(G * (1 - epsilon_L * epsilon_L_F) / (Q + epsilon_L_F * Q_L)) * (1 - np.exp(-(Q + epsilon_L_F * Q_L) * t / V)) if t <= t_g else np.nan for t in time_vector]
    C0 = (G * (1 - epsilon_L * epsilon_L_F)) / (Q + epsilon_L_F * Q_L) * (1 - np.exp(-(Q + epsilon_L_F * Q_L) * t_g / V))
    C_decay = [C0 * np.exp(-(Q + epsilon_L_F * Q_L) * (t - t_g) / V) if t > t_g else np.nan for t in time_vector]
    return time_vector, np.array([x if not np.isnan(x) else y for x, y in zip(C_rise, C_decay)])

def model_110(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, gamma):
    C_avg = (gamma * G * (1 - epsilon_L * epsilon_L_F)) / ((Q + epsilon_RF * Q_R) + epsilon_L_F * Q_L)
    C_LE = C_avg + (epsilon_L * gamma * G) / Q_L
    C_LF = C_LE * (1 - epsilon_L_F)
    C_RF = C_avg * (1 - epsilon_RF)
    return C_avg, C_LE, C_LF, C_RF

def model_111(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, V, t_g, T):
    Q_instead = Q + epsilon_RF * Q_R
    time_vector, concentration = model_109(G, Q_instead, Q_L, epsilon_L, epsilon_L_F, V, t_g, T)
    return time_vector, concentration

# Set the parameters
T = 60
t_g = 15
G = 100
Q = 20
Q_L = 5
epsilon_L = 0.5
epsilon_L_F = 0.75
Q_R = 5
epsilon_RF = 0.9
V = 100
gamma = 0.25  

# Perform calculations
time_109, concentration_109 = model_109(G, Q, Q_L, epsilon_L, epsilon_L_F, V, t_g, T)
time_111, concentration_111 = model_111(G, Q, Q_L, epsilon_L, epsilon_L_F, Q_R, epsilon_RF, V, t_g, T)

# Visualization: Model 109 vs Model 111
plt.figure(figsize=(10, 5))
plt.plot(time_109, concentration_109, label='Model 109')
plt.plot(time_111, concentration_111, label='Model 111', linestyle='--')
plt.title('Model 109 vs Model 111')
plt.xlabel('Time (minutes)')
plt.ylabel('Concentration (mg/m^3)')
plt.legend()
plt.show()

# Calculate Mean Squared Error and R^2
mse = np.mean((concentration_109 - concentration_111)**2)
r2 = np.corrcoef(concentration_109, concentration_111)[0, 1]**2

# Statistical Analysis
normality_test_109 = shapiro(concentration_109[~np.isnan(concentration_109)])
normality_test_111 = shapiro(concentration_111[~np.isnan(concentration_111)])

wilcoxon_test = None
if normality_test_109.pvalue < 0.05 and normality_test_111.pvalue < 0.05:
    wilcoxon_test = wilcoxon(concentration_109[~np.isnan(concentration_109)], concentration_111[~np.isnan(concentration_111)])
    
mse, r2, normality_test_109, normality_test_111, wilcoxon_test

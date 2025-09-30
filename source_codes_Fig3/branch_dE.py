import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

k_agg1 = -45.0
Length = 20e-6
Rc = 0.5e-6
diffL = 20e-9
V_cyt = np.pi * Rc**2 * Length
Na = 6.0223e23 

def func_fit(x, a, b, c):
    #return a * np.exp(x) + b
    return a * x**2 + b * x + c

df = pd.read_csv("./mu-4.5_data.csv")

df = df[df["k_agg1"] == k_agg1]

phi_entire_list = list(np.sort(df["phi_entire"]))

dE_list = []

for pe in phi_entire_list:
    dE_list.append(list(df[df["phi_entire"] == pe]["dE"])[0])

df_extract = pd.DataFrame()
df_extract["conc"] = 1e3 * np.asarray(phi_entire_list) / (Na * V_cyt)
df_extract["phi_entire"] = phi_entire_list
df_extract["dE"] = dE_list

df_extract.to_csv("./K" + str(k_agg1) + "_dE.csv")
popt, pcov = curve_fit(func_fit, df_extract["conc"], df_extract["dE"])
fa = popt[0]
fb = popt[1]
fc = popt[2]
print("Fit values: ", fa, fb, fc)
dE_fit = func_fit(df_extract["conc"], fa, fb, fc)
print(dE_fit)
plt.plot(df_extract["conc"], dE_fit)
plt.plot(df_extract["conc"], df_extract["dE"])
plt.show()

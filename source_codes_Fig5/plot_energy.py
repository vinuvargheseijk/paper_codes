import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import seaborn as sns
sns.set_context("poster")



def func(x, a, b):
    return a * x**2 + b

Length = 10e-6
dendRad = 0.5e-6
V_cyt = np.pi * dendRad**2 * Length
Na = 6.0223e23

mu_list = [-4.0]
phi_entire_list = [200.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0]
agg_list = [-40.0]

for mu in mu_list:
    df = pd.read_csv("./mu" + str(mu) + "_data.csv") 
    #print(df)
    for agg in agg_list:
      energyVsagg = []  
      for phi_entire in phi_entire_list:
         for i in range(len(df)):
            if df.iloc[i]["k_agg1"] == agg:
                if df.iloc[i]["phi_entire"] == phi_entire:
                    print(agg, phi_entire, df.iloc[i]["dE"])
                    energyVsagg.append(df.iloc[i]["dE"])
      conc_list = np.asarray(phi_entire_list) / (V_cyt * Na)      
      #plt.plot(conc_list, energyVsagg, label = "k_agg: " + str(agg))
popt, pcov = curve_fit(func, conc_list, energyVsagg)
print("Curve fit: ", popt)
params = pd.DataFrame()
a = popt[0]
b = popt[1]
print(a,b)
params["p"] = [a,b]
params.to_csv("./params" + str(agg_list[0]) + ".csv")
"""
plt.plot(conc_list, func(conc_list, popt[0], popt[1]), label = "Fit")
plt.title("$\mu_{0}$ = " +str(mu))
plt.xlabel("$\phi_{entire}$")
plt.ylabel("dE")
plt.legend()
plt.show()  
"""
                                             



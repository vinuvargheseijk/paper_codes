import extract_Emin
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
sns.set_context("poster")

phi_entire_list = [2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0]
k_agg = -65.0
wave_sol_dir = "wave_sol_trial"
ws_energy_list = []
emin_energy_list = []
for pe in phi_entire_list:
       file = "mu-4.5_probability_" + str(k_agg) + "_" + str(pe) + "_5.0nm.csv"
       emin, phitot_min = extract_Emin.extract_minEm(pe, file)
       emin_energy_list.append(emin)
       print(pe)
       Ws_file = wave_sol_dir + "/Ws_" + str(pe) + ".csv"
       df_wave = pd.read_csv(Ws_file)
       phitot_ws = float(df_wave[df_wave["Agg. coeff."] == k_agg]["Phitot"])
       energy_ws = float(df_wave[df_wave["Agg. coeff."] == k_agg]["Energy"])
       ws_energy_list.append(energy_ws)
       print("Emin: ", emin, phitot_min)
       print("EWs: ", energy_ws, phitot_ws)
plt.plot(phi_entire_list, np.asarray(ws_energy_list) / 1e3, marker = "*", label = "Ws")
plt.plot(phi_entire_list, np.asarray(emin_energy_list) / 1e3, marker = "o", label = "E_{min}")
plt.xlabel("$\phi_{entire}$")
plt.ylabel("Energy ($10^{3} K_{B}T$)")
plt.legend()
plt.show()


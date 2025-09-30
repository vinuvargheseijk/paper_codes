import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Length = 20e-6
Rc = 0.5e-6
Volume = np.pi * Rc**2 * Length
Na = 6.0223e23
select_kagg = -60.0

phi_entire_list = [2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0]
def plot_WaveVsminEn():
   minEnergy_conc = []
   minEnergy_energy = []
   minEnergy_phitot = []
   for pe in phi_entire_list:
     directory = "../"
     conc = 1e3 * pe / (Na * Volume)
     df = pd.read_csv("./" + directory + "mu-4.5_probability_" + str(select_kagg) + "_" + str(pe) + ".csv")
     min_phitot = list(df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["phitot"])[0]
     min_energy = list(df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["highEnergy"])[0]
     minEnergy_conc.append(conc)
     minEnergy_energy.append(min_energy)
     minEnergy_phitot.append(min_phitot)
   df_calc = pd.read_csv("../wave_solution_final_multi_opt_L20/minData.csv")
   print(minEnergy_phitot)
   return [minEnergy_conc, minEnergy_energy], [df_calc["conc"], df_calc["energy"]]

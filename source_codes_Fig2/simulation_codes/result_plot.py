import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tag
import totalEnergy as TE
import seaborn as sns

sns.set_context("poster")

phi_entire = tag.phi_entire
mem_diffConst = tag.mem_diffConst
mu0_exp = tag.mu0_exp
spacing = tag.spacing
cyt_diffConst = tag.cyt_diffConst
k_agg = tag.k_agg

def_values = [ ["phi_entire", phi_entire], ["num_spines", 1], ["mem_diffConst", mem_diffConst], ["mu0", mu0_exp],["spacing", spacing], ["cyt_diffConst", cyt_diffConst], ["k_agg", k_agg]]

params_list = {"phi_entire":1, "num_spines":1, "mem_diffConst":1, "mu0":1, "spacing":1, "cyt_diffConst":1, "k_agg":1}
def extract_data(params_vary, value_range, N, phi_entire = phi_entire):
  ns = N
  params_list_const = params_list.copy() 
  print(params_list_const)
  params_list_const.pop(params_vary)  
  print("Constant: ", params_list_const)
  p_count = 0
  for p in params_list:
      if p in params_list_const:
          params_list[p] = def_values[p_count][1]
          print(params_list[p])
      else:
          params_list[p] = value_range
          print(p)
      p_count = p_count + 1    
  print("After modifying: ", params_list)    
  #tag_string = "phiE_" + str(params_list["phi_entire"]) + "_ns_" + str(params_list["num_spines"]) + "_mD_" + str(params_list["mem_diffConst"]) + "_mu0_" + str(params_list["mu0"]) + "_spacing_" + str(params_list["spacing"]) + "_cD_" + str(params_list["cyt_diffConst"]) + "_k_agg_" + str(params_list["k_agg"])
  tag_string = "phiE_" + str(phi_entire) + "_ns_" + str(params_list["num_spines"]) + "_mD_" + str(params_list["mem_diffConst"]) + "_mu0_" + str(params_list["mu0"]) + "_spacing_" + str(params_list["spacing"]) + "_cD_" + str(params_list["cyt_diffConst"]) + "_k_agg_" + str(params_list["k_agg"])
  print(tag_string)
  tEnergy = 0
  for n in range(1, ns+1):  
    phi_df = pd.read_csv(tag_string + "phi.csv")
    hmean = pd.read_csv(tag_string + "hmean.csv")
    rp_df = pd.read_csv(tag_string + "rp.csv")
    theta_df = pd.read_csv(tag_string + "theta.csv")
    conc_df = pd.read_csv(tag_string + "cConc.csv")
    rm = 1 / hmean.iloc[-1][n]
    theta = theta_df.iloc[-1][n]
    rp = rp_df.iloc[-1][n]
    r0 = (rm + rp) * np.sin(theta)
    conc = conc_df.iloc[-1][n]
    if np.isnan(theta) == False:
      area = 2 * np.pi * rm * rm * (1 - np.cos(theta))
      phitot = phi_df.iloc[-1][n] * area * TE.phisat
      current_energy = TE.total_energy([rm * 1e6, theta, rp * 1e6], phitot, TE.phisat, TE.khat, TE.fs, TE.k, 1, 0, 0, 0, conc) * 1e-16 / tag.KBT
      tEnergy = tEnergy + current_energy
      print("Current energy: ", current_energy)
      print("Current phitot: ", phitot)
  print("Total Energy of " + str(ns) + " spines: ", tEnergy)  
  return np.asarray(r0) * 1e6, tEnergy



fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3,2, figsize = (12,12))
phi_entire_range = [4000, 6000]
mu0_range = [-4.0, -4.5, -5.0]
for phiE in phi_entire_range:
  r0 = []
  for i in mu0_range:
    r0i, tEn = extract_data("mu0", i, 1, phiE)  
    r0.append(r0i)
  ax1.set_ylabel("$r_{0} (\mu m)$")
  ax1.set_xlabel("$\mu_{0}$")
  ax1.plot(mu0_range, r0, marker = "o", label = "$\phi_{entire}$ = " + str(phiE))  #list reversed to show increasing attractive forces
ax1.invert_xaxis()
ax1.legend()
ax1.set_yticks([0.0, 0.5, 1.0])
input()

n_range = [1, 2, 3, 4, 5, 6]
for phiE in phi_entire_range:
  r0 = []
  tEn_list = []
  for i in n_range:
    r0i, tEn = extract_data("num_spines", i, i, phiE)  
    r0.append(r0i)
    tEn_list.append(tEn)
  ax2.plot(n_range, r0, marker = "o")
  ax2.set_ylabel("$r_{0} (\mu m)$")
  ax2.set_xlabel("$N$")
  ax2.set_yticks([0.0, 0.5, 1.0])
  ax2.set_title("$\mu_{0}$ = " + str(mu0_exp))

  ax3.plot(n_range, tEn_list, marker = "o")
  ax3.set_ylabel("Energy")
  ax3.set_xlabel("$N$")
  ax3.set_title("$\mu_{0}$ = " + str(mu0_exp))
              
"""
E_fr_data = pd.read_csv("./energy_fraction.csv")
fraction_data = []
for i in range(1, len(E_fr_data.iloc[0]) - 1):
    fraction_data.append(E_fr_data.iloc[0][i])
print(len(E_fr_data.columns.values[1:-1]))
print(fraction_data)
"""
#ax4.pie(np.asarray(fraction_data), labels = E_fr_data.columns.values[1:-1])
plt.show()

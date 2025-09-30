import pandas as pd
import numpy as np
import writeXML
import matplotlib.pyplot as plt
import readXML as RX
import xml.etree.ElementTree as ET
import seaborn as sns
from scipy import optimize
import energy_process as EP
import argparse
sns.set_context("poster")

parser = argparse.ArgumentParser()
parser.add_argument("--mu0", type = float)
args = parser.parse_args()

KBT = EP.KBT
Na = EP.Na
RT = EP.RT
RT_kcal_mol = EP.RT_kcal_mol #kcal/mol
k = EP.k
khat = EP.khat
kentropy = KBT
Cp = EP.Cp
phisat = EP.phisat
phi_entire_list = []
k_agg_list = []
state_list = []
sum_list = []
energy_list = []
dE_list = []
phi_list = []
diff_phitot_list = []
phitot_list = []
prot_size = []
Cp1_list = []
Cp2_list = []
scalF = 1.0
Kb = 40.0
spacing = 2.0e-6
dt = 0.24
shape_dt = 0.2
plot_time = 299
Length = 20e-6
Rc = 0.5e-6
Na = 6.023e23
V_cyt = np.pi * Rc**2 * Length 

fig, (ax1, ax2) = plt.subplots(2)
phi_entire_list = [2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0]
agg_list = [-45.0, -50.0, -55.0, -60.0, -65.0]
mu0_list = [-4.5]
spacing_list = np.asarray([0.75, 1.0, 1.5, 3.0, 5.0])
scaling_list = [0.15, 0.25, 0.35, 0.45, 0.55, 0.75, 0.85]
n_list = [1]
cyt_diff_list = np.asarray([0.4]) * 1e-12
cyt_diffConst = cyt_diff_list[0]
var_param1 = phi_entire_list
const_param = agg_list[0]
const_param_name = "Aggregation"
var_param2 = n_list
const_param2 = mu0_list[0]
const_param2_name = "mu0"
direction = "forward"

hmp_r1 = ["phi_entire"] #Check if const and var params are correct
hmp_r2 = ["N"]

def line_plots(v_parameters, c_parameters, df, plot_name):
    print(df[c_parameters[0]][0])
    plot_list = []
    for i in range(len(df)):
        print(df[c_parameters[0]][i], c_parameters[1])
        if df[c_parameters[0]][i] == c_parameters[1]:
            print("Found: ", df[c_parameters[0]][i])
            for j in v_parameters[1]:
               print(j) 
               if df[v_parameters[0]][i] == j:
                   print("Found", df[v_parameters[0]][i])
                   print(df[plot_name][i])
                   plot_list.append(df[plot_name][i])
    return plot_list            

def split_list(list_name):
    len_split = len(hmp_r1)
    split_list = []
    count = 0
    for l in range(len(var_param1)):
        split_list.append(list_name[count : count + len(var_param2)])
        count = count + len(var_param2)
    print(split_list)   
    return split_list



for vp1 in var_param1:
    for vp2 in var_param2:
       tag_string = "phi_entire_" + str(vp1) + "_k_agg_" + str(const_param) + "_ns_" + str(vp2) + "_mu0_" + str(const_param2) + "_spacing_" + str(spacing) + "_dt_" + str(dt) + "_cytD_" + str(cyt_diffConst)
       #in_index = -1
       in_index = int(plot_time / shape_dt)
       print("TAG STRING: ", tag_string)
       hmp_r1.append(vp1)
       hmp_r2.append(vp2)
       Cp1_list.append(const_param)
       Cp2_list.append(const_param2)

       try:
          print("reading file") 
          phi = pd.read_csv(tag_string + "phi.csv")
          theta = pd.read_csv(tag_string + "theta.csv")
          saddle_start = pd.read_csv(tag_string + "sstart.csv")
          hmean = pd.read_csv(tag_string + "hmean.csv")
          rp = pd.read_csv(tag_string + "rp.csv")
          mConc = pd.read_csv(tag_string + "mConc.csv")
          cConc = pd.read_csv(tag_string + "cConc.csv")
          total_energy = 0
          total_phitot = 0
          state = ""
          sum_t = 0
          phitot_ns = []
          #counting spines
          sp_count = 0
          for ns in range(1, 10):
              try:
                  phi.iloc[in_index][ns]
                  sp_count = sp_count + 1
              except:
                  print("Num spines: ", sp_count)
          for ns in range(1, sp_count+1):
             try:  
               conc = cConc.iloc[in_index][ns]
               rp_last = rp.iloc[in_index][ns]
               rm_last = 1 / hmean.iloc[in_index][ns]
               theta_last = theta.iloc[in_index][ns]
               init_area = EP.calc_dome_area(1/hmean.iloc[in_index][ns], theta.iloc[in_index][ns])
               area = EP.calc_dome_area(rm_last, theta_last)
               phi_last = phi.iloc[in_index][ns]
               init_phitot = phi.iloc[in_index][ns] * phisat * init_area
               phitot = phi_last * phisat * area
               prot_size_i = 2 * (rp_last + rm_last) * np.sin(theta_last)
               phitot_ns.append(phitot)
               init_energy = EP.total_energy([(1/hmean.iloc[in_index][ns]) * 1e6 * theta.iloc[in_index][ns], theta.iloc[in_index][ns], rp.iloc[in_index][ns] * 1e6], init_phitot, cConc.iloc[in_index][ns], const_param, const_param2)  #Check for input arguments
               energy = EP.total_energy([rm_last * 1e6 * theta_last, theta_last, rp.iloc[in_index][ns] * 1e6], phitot, conc, const_param, const_param2) #Check for input arguments
               if np.isnan(energy) == False:
                 total_energy = total_energy + energy
                 total_phitot = total_phitot + phitot
                 
               if phi.iloc[in_index][ns] <= 0.1:
                   state = state + "0"
                   sum_t = sum_t + 0
               else:
                   state = state + "1"
                   sum_t = sum_t + 1
             except:
                   print("Number of spines: " + str(ns))
          energy_list.append(total_energy)
          state_list.append(state)
          phitot_list.append(total_phitot)
          prot_size.append(prot_size_i)
          try:
             diff_phitot_list.append(phitot_ns[0] - phitot_ns[1])
          except:
             diff_phitot_list.append(np.nan)  
          sum_list.append(sum_t)
          phi_list.append(phi.iloc[in_index][1])
       except Exception as excpt:
           type_excpt = type(excpt).__name__
           print(type_excpt)
           print(tag_string)
           if type_excpt != "IndexError":
             energy_list.append(np.nan)
             phitot_list.append(np.nan)
             state_list.append(np.nan)
             diff_phitot_list.append(np.nan)
             phi_list.append(np.nan)
             sum_list.append(np.nan)
             prot_size.append(np.nan)
               

split_list(state_list)
df = pd.DataFrame()
print("Energy: ", energy_list)
print("phitot diff: ", diff_phitot_list)
print(len(energy_list))
print(len(diff_phitot_list))
print(len(phi_entire_list))
df[hmp_r1[0]] = hmp_r1[1:]
df[hmp_r2[0]] = hmp_r2[1:]
df["state"] = state_list
df["energy"] = energy_list
df["phi"] = phi_list
df["diff_phitot"] = diff_phitot_list
df["num_final"] = sum_list
df["phitot"] = phitot_list
df["protSize"] = prot_size
df[const_param_name] = Cp1_list
df[const_param2_name] = Cp2_list
print(state_list)
print(df)


print("Data")
data = np.array(split_list(diff_phitot_list))
print("Labels")
labels = np.array(split_list(state_list))
res = df.pivot(index = hmp_r1[0], columns = hmp_r2[0], values = "phitot")
ax1.set_title("$\phi_{tot}$, " + const_param_name + str(const_param))
ax1 = sns.heatmap(res, cmap = "cividis", ax = ax1, annot = True)
ax1.invert_xaxis()
ax1.invert_yaxis()
"""
res = df.pivot(index = hmp_r1[0], columns = hmp_r2[0], values = "diff_phitot")
ax1.set_title("$\phi_{tot1} - \phi_{tot2}$, " + const_param_name + str(const_param))
ax1 = sns.heatmap(res, cmap = "cividis", ax = ax1, annot = True)
ax1.invert_xaxis()
ax1.invert_yaxis()


ax1 = sns.heatmap(data, ax= ax1, annot = labels, xticklabels = scaling_list, yticklabels = phi_entire_list, fmt = '')
ax1.invert_xaxis()
ax1.invert_yaxis()
"""
res_e = df.pivot(index = hmp_r1[0], columns = hmp_r2[0], values = "energy")
ax2.set_title("Energy, " + const_param_name + str(const_param))
ax2 = sns.heatmap(res_e, cmap = "cividis", ax = ax2, annot = True, fmt = ".3f")
ax2.invert_xaxis()
ax2.invert_yaxis()

try:
   pd.read_csv("./sim.csv") 
   df.to_csv("./sim.csv", mode = 'a', index = False, header = False)
except:   
   df.to_csv("./sim.csv", index = False)
fig, (ax1, ax2) = plt.subplots(2)
for pe in phi_entire_list:
  cyt_conc = (pe / (V_cyt * Na)) * 1e3  
  plot_list = line_plots([hmp_r2[0], n_list], [hmp_r1[0], pe], df, "energy")
  print("PLOT LIST: ", plot_list)
  ax1.plot(n_list, -1 * np.asarray(plot_list) / min(plot_list), marker = "*", label = "$Conc. (\mu M)$" + str(round(cyt_conc,1)))
ax1.legend()  
ax1.set_xlabel("N protrusions")
ax1.set_ylabel("Norm. energy")
ax1.set_title("Diffusion $\mu m^{2}/s$  = " + str(cyt_diffConst * 1e12))

"""
for pe in phi_entire_list:
  plot_list = line_plots([hmp_r2[0], spacing_list], [hmp_r1[0], pe], df, "diff_phitot")
  print("PLOT LIST: ", plot_list)
  ax2.plot(spacing_list, np.asarray(plot_list), marker = "*", label = "$\phi_{entire}$" + str(pe))
ax2.legend()  :q

ax2.set_xlabel("Spacing ($\mu m$)")
ax2.set_ylabel("diff. $\phi_{tot}$")
"""

plt.show()


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import readXML as RX
import xml.etree.ElementTree as ET
import seaborn as sns
import energy_process as EP
from mpl_toolkits.mplot3d import Axes3D
sns.set_context("poster")

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
global state_list, energy_list, dE_list, phi_list, diff_phitot_list, phitot_list
state_list = []
sum_list = []
energy_list = []
dE_list = []
phi_list = []
diff_phitot_list = []
phitot_list = []
Ca_list = []
actCaMKII_list = []
Tiam_list = []
Rac_list = []
cytosol_list = []
mu0_exp = -4.5
scalF = 1.0
Kb = 40.0
dt = 0.24
cyt_diffConst = 0.5e-12
Length = 10e-6
diffL = 20e-9
Rc = 0.5e-6
Na = 6.023e23
V_cyt = np.pi * Rc**2 * Length 
V_voxel = np.pi * 0.5 * (1e-6)**2 * diffL   
shape_dt = 0.2
CaTau = 0.08
stimTypes = ["s", "b"]
num_bursts = {"s":30, "b":3}
#freq_list = [0.5, 1.0, 2.0]
freq_list = [0.5]
ypos_list = [0.1e-6, 0.25e-6, 0.5e-6, 0.75e-6, 1.5e-6, 2.0e-6]
source = 8000.0
num_spines = 1
y_pos = 0.25e-6
spacing = 0.75e-6

fig1 = plt.figure(figsize=(12,10))
ax00 = fig1.add_subplot(121, projection='3d')
ax00.view_init(elev=-30, azim=45, roll=90)
ax01 = fig1.add_subplot(122, projection='3d')
ax01.view_init(elev=-30, azim=45, roll=90)
fig2 = plt.figure(figsize=(12,10))
ax1 = fig2.add_subplot(121, projection='3d')
ax2 = fig2.add_subplot(122, projection='3d')
ax1.view_init(elev=-30, azim=45, roll=90)
ax2.view_init(elev=-30, azim=45, roll=90)
fig3 = plt.figure(figsize=(12,10))
ax3 = fig3.add_subplot(121, projection='3d')
ax4 = fig3.add_subplot(122, projection='3d')
ax3.view_init(elev=-30, azim=45, roll=90)
ax4.view_init(elev=-30, azim=45, roll=90)
fig4 = plt.figure(figsize=(12,10))
ax5 = fig4.add_subplot(121, projection='3d')
ax6 = fig4.add_subplot(122, projection='3d')
ax5.view_init(elev=-30, azim=45, roll=90)
ax6.view_init(elev=-30, azim=45, roll=90)
fig5 = plt.figure(figsize=(12,10))
ax7 = fig5.add_subplot(121, projection='3d')
ax8 = fig5.add_subplot(122, projection='3d')
ax7.view_init(elev=-30, azim=45, roll=90)
ax8.view_init(elev=-30, azim=45, roll=90)


settling_time = 15.0
for freq in range(len(freq_list)):
  burst_plot_axes = [ax01, ax2, ax4, ax6, ax8]
  single_plot_axes = [ax00, ax1, ax3, ax5, ax7]
  for stimT in stimTypes:
       tag_string = "Source_" + str(source) + "_nb_" + str(num_bursts[stimT]) + "_ns_" + str(num_spines) + "_ypos_" + str(y_pos) + "_spacing_" + str(spacing) + "_freq_" + str(freq_list[freq]) + "_CaTau_" + str(CaTau) + "_typeStim_" + str(stimT)
       chemTimes = pd.read_csv(tag_string + "chemTimes.csv")
       settling_index = np.where(chemTimes["time"] > 15.0)[0][0] - 1
       print(settling_index)
       chemTimes = chemTimes["time"]
       chemTimes_settle = chemTimes
       CalciumXML = tag_string + "Ca.xml"
       actCaMKIIXML = tag_string + "thr286.xml" 
       TiamXML = tag_string + "Tiam.xml"
       RacXML = tag_string + "Rac.xml"
       cytosolXML = tag_string + "cytosol.xml"

       check_files = [CalciumXML, actCaMKIIXML, TiamXML, RacXML, cytosolXML]
       mol_names = ["Calcium", "$CaMKII_{act}$", "Tiam", "Rac", "$IRSp53_{act}$"] 
       mol_count = 0
       for c in check_files:
          filename = c
          if filename == CalciumXML:
            plot3_color = "m"
          if filename == cytosolXML:
            plot3_color = "g"
          else:
            plot3_color = "b"
          print(filename)
          tree1 = ET.parse(filename)
          tree = [tree1]
          sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
          time_axis = list(np.asarray(chemTimes.iloc[settling_index:-1]) - settling_time)
          plot_times = []
          plot_indices = []
          sum_conc_indices = []
          print(time_axis)
          print(len(time_axis))
          conc_mols = 1e3 * np.asarray(sum_an) / (V_cyt * Na)
          conc_mols = conc_mols[settling_index:-1]
          spatial_mols = an_y[settling_index:-1]
          x_axis = np.linspace(0, Length * 1e6, len(spatial_mols[0]))
          #Finding max value in all an_y
          max_number = max(sum_an)
          max_conc = 1e3 * max_number / (V_cyt * Na)
          if stimT == "s": 
             plot3_interval = 10
             for ta in range(0, len(time_axis), plot3_interval):
                 plot_times.append(time_axis[ta])
                 plot_indices.append(ta)
                 sum_conc_indices.append(1e3 * sum_an[ta + settling_index] / (V_cyt * Na))
             for p3t in range(len(plot_times)):
                  spatial_conc = 1e3 * np.asarray(spatial_mols[plot_indices[p3t]]) / (V_voxel * Na)
                  if sum_conc_indices[p3t] / max_conc > 1.0:
                     input()
                  single_plot_axes[0].plot3D(x_axis, spatial_conc, time_axis[plot_indices[p3t]], plot3_color, alpha = sum_conc_indices[p3t] / max_conc )
                  single_plot_axes[0].set_xlabel("x $\mu m$")
                  single_plot_axes[0].set_ylabel("$\mu M$")
                  single_plot_axes[0].set_zlabel("Time (s)")
                  single_plot_axes[0].set_title(mol_names[mol_count])
                  single_plot_axes[0].set_zticks([0, 100, 200])
             try:
                single_plot_axes.pop(0)
             except:
                print("Plotting of " + str(stimT) + " ended")
          if stimT == "b":
             if filename == CalciumXML:
                 plot3_interval = 1
             else:
                 plot3_interval = 10
             for ta in range(0, len(time_axis), plot3_interval):
                 plot_times.append(time_axis[ta])
                 plot_indices.append(ta)
                 sum_conc_indices.append(1e3 * sum_an[ta + settling_index] / (V_cyt * Na))
             for p3t in range(len(plot_times)):
                  spatial_conc = 1e3 * np.asarray(spatial_mols[plot_indices[p3t]]) / (V_voxel * Na)
                  if sum_conc_indices[p3t] / max_conc > 1.0:
                     input()
                  burst_plot_axes[0].plot3D(x_axis, spatial_conc, time_axis[plot_indices[p3t]], plot3_color, alpha = sum_conc_indices[p3t] / max_conc )
                  burst_plot_axes[0].set_xlabel("x $\mu m$")
                  burst_plot_axes[0].set_zticks([0, 100, 200])
                  #burst_plot_axes[0].set_ylabel("$\mu M$")
                  #burst_plot_axes[0].set_zlabel("Time (s)")
             try:
                burst_plot_axes.pop(0)
             except:
                print("Plotting of " + str(stimT) + " ended")
          mol_count = mol_count + 1

fig2, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2)
ypos_df = []
csmax_df = []
cbmax_df = []
rsmax_df = []
rbmax_df = []
for ypos in ypos_list:
  print(ypos)
  burst_plot_axes = [ax2, ax4]
  single_plot_axes = [ax1, ax3]
  for stimT in stimTypes:
       tag_string = "Source_" + str(source) + "_nb_" + str(num_bursts[stimT]) + "_ns_" + str(num_spines) + "_ypos_" + str(ypos) + "_spacing_" + str(spacing) + "_freq_" + str(0.5) + "_CaTau_" + str(CaTau) + "_typeStim_" + str(stimT)
       chemTimes = pd.read_csv(tag_string + "chemTimes.csv")
       settling_index = np.where(chemTimes["time"] > 15.0)[0][0] - 1
       print(settling_index)
       chemTimes = chemTimes["time"]
       chemTimes_settle = chemTimes
       CalciumXML = tag_string + "Ca.xml"
       actCaMKIIXML = tag_string + "thr286.xml" 
       TiamXML = tag_string + "Tiam.xml"
       RacXML = tag_string + "Rac.xml"
       cytosolXML = tag_string + "cytosol.xml"
       filenames = [CalciumXML, cytosolXML]
       for filename in filenames:
         tree1 = ET.parse(filename)
         tree = [tree1]
         sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
         conc_mols = 1e3 * np.asarray(sum_an) / (V_cyt * Na)
         conc_mols = conc_mols[settling_index:-1]
         time_axis = list(np.asarray(chemTimes.iloc[settling_index:-1]) - settling_time)
         if stimT == "s":
            if filename == CalciumXML:
               ypos_df.append(ypos * 1e6)
               csmax_df.append(1e3 * max(sum_an) / (V_cyt * Na))
            if filename == cytosolXML:
               rsmax_df.append(1e3 * max(sum_an) / (V_cyt * Na))
         if stimT == "b":
            if filename == CalciumXML:
               cbmax_df.append(1e3 * max(sum_an) / (V_cyt * Na))
            if filename == cytosolXML:
               rbmax_df.append(1e3 * max(sum_an) / (V_cyt * Na))

ax1.set_ylabel("$\mu M$")
ax2.plot(ypos_df, cbmax_df, "m", marker = "*")
ax3.plot(ypos_df, rsmax_df, "g", marker = "*")
ax3.set_ylabel("$\mu M$")
ax3.set_xlabel("y $\mu m$")
ax4.plot(ypos_df, rbmax_df, "g", marker = "*")
ax4.set_xlabel("y $\mu m$")


plt.show()


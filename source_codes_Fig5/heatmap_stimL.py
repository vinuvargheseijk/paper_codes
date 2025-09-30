import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import readXML as RX
import xml.etree.ElementTree as ET
import seaborn as sns
import energy_process as EP
import plot_shape
from mpl_toolkits.mplot3d import Axes3D
sns.set_context("paper")


KBT = EP.KBT
Na = EP.Na
RT = EP.RT
RT_kcal_mol = EP.RT_kcal_mol #kcal/mol
k = EP.k
khat = EP.khat
kentropy = KBT
Cp = EP.Cp
typeStim = "b"
yratio = 1.0
ypos = 0.25
nb_perturb = 4.0
num_spines = 9
spacing = 1.0
stimL_list = [6.0, 8.0, 10.0]
freq_list = [1.0]
nb = 8

def plot_axis(x_length, y_length):
    axis = ax2.inset_axes([-0.1, 0, 0.15, 0.15])
    axis.plot([0, x_length], [0, 0 * 1e-3], 'k')
    axis.plot([0, 0], [0, y_length], 'k')
    axis.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlim([0, 20])
    axis.set_ylim([0, 55])
    #axis.set_xlabel(str(x_length) + "$\mu m$")
    axis.set_xlabel(str(x_length) + "$\mu m$", loc="left", rotation="horizontal")
    axis.set_ylabel(str(y_length) + "nm")

def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

fig = plt.figure(figsize = (10, 12))
ax1 = plt.subplot2grid((2, 2), (0, 0), colspan = 1) 
ax11 = plt.subplot2grid((2, 2), (0, 1), colspan = 1) 
ax2 = plt.subplot2grid((2, 2), (1, 0), colspan = 2) 

inset_count_y = 0.8
init_flag = True
for stimL in stimL_list:
    num_list = []
    cluster_width = []
    inset_count_x = 0.1
    globals()["stimL NS" + str(stimL)] = []
    globals()["stimL CW" + str(stimL)] = []
    for fq in freq_list:
       tag_string = "typeStim_" + str(typeStim) + "_nb_" + str(nb) + "_ns_" + str(num_spines) + "_stimL_" + str(stimL) + "_spacing_" + str(spacing) + "_freq_" + str(fq) + "_ypos_" + str(ypos) + "_nb_perturb_" + str(nb_perturb)
       print(tag_string)
       phi = pd.read_csv("./" + tag_string + "phi.csv")
       hmean = pd.read_csv("./" + tag_string + "hmean.csv")
       theta = pd.read_csv("./" + tag_string + "theta.csv")
       saddle_start = pd.read_csv(tag_string + "sstart.csv")
       rp = pd.read_csv(tag_string + "rp.csv")
       try:
          last_phi = remove_nan(phi.iloc[-1][1:])
       except:
          globals()["stimL NS" + str(stimL)].append(np.nan)  
          globals()["stimL CW" + str(stimL)].append(np.nan)
          num_list.append(np.nan) 
          cluster_width.append(np.nan)
          continue 
       last_hmean = remove_nan(hmean.iloc[-1][1:])
       last_rp = remove_nan(rp.iloc[-1][1:])
       last_theta = remove_nan(theta.iloc[-1][1:])
       last_saddle_start = remove_nan(saddle_start.iloc[-1][1:])
       print(last_phi)
       spine_number = len(last_phi)
       globals()["stimL NS" + str(stimL)].append(spine_number)  
       globals()["stimL CW" + str(stimL)].append(last_saddle_start[-1] - last_saddle_start[0])
       num_list.append(spine_number)
       cluster_width.append(last_saddle_start[-1] - last_saddle_start[0])
       if init_flag and nb == 10 and fq == 0.5:
         axin_init = ax1.inset_axes([0.7, 0.3, 0.15, 0.15])
         pos = 90
         print(saddle_start.iloc[pos][1:])
         init_x, init_y = plot_shape.get_shape(remove_nan(saddle_start.iloc[pos][1:]), remove_nan(theta.iloc[pos][1:]), remove_nan(hmean.iloc[pos][1:]), remove_nan(rp.iloc[pos][1:]))
         axin_init.plot(init_x, np.asarray(init_y) * 1e3)
         init_flag = False
       if len(last_phi) > 0:
          x, y = plot_shape.get_shape(last_saddle_start, last_theta, last_hmean, last_rp)
          axin1 = ax2.inset_axes([inset_count_x, inset_count_y, 0.15, 0.15])
          inset_count_x = inset_count_x + 0.3
          axin1.plot(x, np.asarray(y) * 1e3)
          axin1.set_ylabel("(" + str(fq) + ", " + str(nb) + ")")
          axin1.set_xlim([0, 20])
          axin1.set_ylim([1.2, 55])
          axin1.set_yticks([1.2, 55])
          plot_axis(2, 50)
          axin1.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
          axin1.set_xticks([])
          axin1.set_yticks([])
          axin1.plot([10] * 5, [0, 5, 10, 20, 40], linestyle = "--", color = 'k')
          axin1.plot([stimL] * 5, [0, 5, 10, 20, 40], linestyle = "-", color = 'b')
          if fq != 4.0:
             axin1.set_xticks([])
       else:   
          axin1 = ax2.inset_axes([inset_count_x, inset_count_y, 0.15, 0.15])
          axin1.plot()
          inset_count_x = inset_count_x + 0.3
          axin1.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
          axin1.set_xticks([])
          axin1.set_yticks([])
       inset_count_y = inset_count_y - 0.4
for stimL in stimL_list:       
  ax1.plot(freq_list, globals()["stimL NS" + str(stimL)], marker = "*", label = "sl: " + str(stimL))   
  ax11.plot(freq_list, np.asarray(globals()["stimL CW" + str(stimL)]) * 1e6, marker = "*", label = "sl: " + str(stimL))   
ax1.set_ylabel("# spines")    
ax11.set_ylabel("Cluster width $\mu m$")    
ax1.set_xlabel("# bursts total")    
ax1.set_yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
ax1.legend() 
ax11.set_ylim([0, 20])
ax11.legend()
ax2.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
ax2.set_xticks([])
ax2.set_yticks([])
plt.show()

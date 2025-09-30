import pandas as pd
import dendShape
import math
import numpy as np
import matplotlib.pyplot as plt
import readXML as RX
import xml.etree.ElementTree as ET
import seaborn as sns
import energy_process as EP
import shankReac as sReac
import plot_shape
from mpl_toolkits.mplot3d import Axes3D
import readXML as RX
import NVsConc_energy as NCE
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle 
from matplotlib.colors import LogNorm, Normalize
from matplotlib.colors import LinearSegmentedColormap
sns.set_context("poster")
import pylustrator
pylustrator.start()



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
ypos_list = [0.25]
nb_perturb = 4.0
num_spine_list = [2, 3]
cyt_conc_list = np.asarray([0.6e-3, 0.8e-3, 1.0e-3, 1.2e-3]) * 1e3
spine_spacing_list = [0.6, 0.7, 0.8, 1.0, 1.5, 2.0, 2.5]
spine_create_delay_list = [2.0, 5.0, 10.0, 15.0]
time_of_plot = 60
trial_plots = True
phisat = EP.phisat
Length = 20e-6
diffL = 20e-9
v_voxel = np.pi * (EP.Rc)**2 * diffL
v_cyt = np.pi * (EP.Rc)**2 * Length
plot_save_dt = 0.05
shape_dt = 0.2
threshold_val = 0.0002
membrane_detect = 1e-1

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, s = 0.5)

    # now determine nice limits by hand:
    binwidth = 0.05
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')

def get_parameters(conc, spacing, spine_index, delay):
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "hmean.csv")
   print("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "hmean.csv")
   print(df)
   for i in range(len(df)):
      if np.isnan(df.iloc[i]["1"]) == False:
         init_loc = i
         break
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "shapeTimes.csv")
   print("Time: ", df.iloc[init_loc]["0"])
   time_2_init = df.iloc[init_loc]["0"]
   params = ["hmean", "theta", "sstart", "rp"]
   for p in params:
      df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + p + ".csv")
      globals()[p] = df.iloc[init_loc][str(spine_index)]
      print(globals()[p])
   x, y = dendShape.get_shape([globals()["sstart"]], [globals()["theta"]], [globals()["hmean"]], [globals()["rp"]])
   return x, y, time_2_init

def get_parameters_merge(conc, spacing, spine_index, delay):
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "hmean.csv")
   print(df)
   for i in range(len(df)):
      if np.isnan(df.iloc[i]["1"]) == False:
         init_loc = i
         break
   for i in range(init_loc, len(df)):
      if np.isnan(df.iloc[i]["1"]) == True:
         merge_loc = i - 1
         break
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "shapeTimes.csv")
   print("Time: ", df.iloc[merge_loc]["0"])
   time_2_merge = df.iloc[merge_loc]["0"]
   params = ["hmean", "theta", "sstart", "rp"]
   for p in params:
      df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + p + ".csv")
      globals()[p] = df.iloc[merge_loc][str(spine_index)]
      print(globals()[p])
   x, y = dendShape.get_shape([globals()["sstart"]], [globals()["theta"]], [globals()["hmean"]], [globals()["rp"]])
   return x, y, time_2_merge


def get_parameters_select_time(conc, spacing, time, spine_index, delay):
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "hmean.csv")
   print("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "hmean.csv")
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "shapeTimes.csv")
   loc_index = int(time / 0.2)
   time_2_init = df.iloc[loc_index]["0"]
   params = ["hmean", "theta", "sstart", "rp"]
   for p in params:
      df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + p + ".csv")
      globals()[p] = df.iloc[loc_index][str(spine_index)]
      print(globals()[p])
   x, y = dendShape.get_shape([globals()["sstart"]], [globals()["theta"]], [globals()["hmean"]], [globals()["rp"]])
   return x, y


def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

def min_space(axis):
    N = 2
    conc_list = []
    delay_list = []
    spacing_list = []
    N_list = []
    df = pd.DataFrame()
    cyt_conc_list = np.asarray([0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3, 1.2e-3, 1.4e-3])
    spine_create_delay_list = [0.0, 5.0, 15.0, 25.0]
    for conc in cyt_conc_list:
       for spine_create_delay in spine_create_delay_list:
           minimum_spacing = spine_spacing_list[-1]
           for spine_spacing in spine_spacing_list:
               tag_string = "CytConc_" + str(conc) + "_ns_" + str(N) + "_spine_create_delay_" + str(spine_create_delay) + "_spine_spacing_" + str(spine_spacing)
               phi = pd.read_csv("./" + tag_string + "phi.csv")
               finalN = len(remove_nan(phi.iloc[-1][1:]))
               if finalN > 1:
                  if spine_spacing < minimum_spacing:
                     minimum_spacing = spine_spacing
                     maxN = finalN
           N_list.append(maxN)           
           spacing_list.append(minimum_spacing)
           conc_list.append(conc * 1e3)
           delay_list.append(spine_create_delay)
    df["Total Conc"] = conc_list
    df["delay"] = delay_list
    df["spacing"] = spacing_list
    df["FinalN"] = N_list
    res = df.pivot(index = "Total Conc", columns = "delay", values = "spacing")
    axis = sns.heatmap(res, cmap = 'cividis', ax = axis, cbar_kws = {"label" : "Minimum spacing $\mu m$"})
    axis.set_yticklabels(axis.get_yticklabels(), rotation = 0)
    axis.invert_yaxis()
    return df



def extract_conc_profile(cytConc, N, spine_create_delay, spine_spacing, filename, center, width, zoom_time):
       tag_string = "CytConc_" + str(round(cytConc * 1e-3, 5)) + "_ns_" + str(N) + "_spine_create_delay_" + str(spine_create_delay) + "_spine_spacing_" + str(spine_spacing)
       print("TAG: ", tag_string)
       filename = tag_string + filename + ".xml"
       cytosol = tag_string + "cytosol.xml"
       chemTimes = pd.read_csv(tag_string + "chemTimes.csv")["time"]
       tree1 = ET.parse(filename)
       tree = [tree1]
       tree2 = ET.parse(cytosol)
       tree_c = [tree2]
       time_heat = []
       x_heat = []
       value_heat = []
       x_axis = np.arange(0, Length * 1e6, diffL * 1e6)
       center_voxel = (center) / (diffL * 1e6)
       width_voxel = (width) / (diffL * 1e6)
       sum_an,max_an_y,an_y = RX.plotXML(filename,tree)
       sum_Can,max_Can_y,Can_y = RX.plotXML(cytosol,tree_c)
       df = pd.DataFrame()
       lowerTime = chemTimes[chemTimes > zoom_time[0]]
       UpperTime = chemTimes[chemTimes < zoom_time[1]]
       lowerTime_index = list(lowerTime.index)[0]
       upperTime_index = list(UpperTime.index)[-1]
       print("INDICES: ", lowerTime_index)
       print("INDICES: ", upperTime_index)
       conc1 = []
       conc2 = []
       cytConc1 = []
       cytConc2 = []
       last_time = 0
       shape_count = 0
       recorded_chemTimes = chemTimes[lowerTime_index: upperTime_index]
       for t in range(lowerTime_index, upperTime_index):
           spine1_voxels = []
           spine2_voxels = []
           sum_mem1 = 0
           sum_mem2 = 0
           sum_cyt1 = 0
           sum_cyt2 = 0
           for voxels in range(int(center_voxel - width_voxel), int(center_voxel)):
               if an_y[t][voxels] > membrane_detect:
                 sum_mem1 = sum_mem1 + an_y[t][voxels]
                 sum_cyt1 = sum_cyt1 + Can_y[t][voxels]
                 spine1_voxels.append(voxels)
           if len(spine1_voxels) > 0:
                temp_conc = (sum_cyt1 + sum_mem1) / (Na * v_voxel * len(spine1_voxels))
                temp_cytConc = (sum_cyt1) / (Na * v_voxel * len(spine1_voxels))
                cytConc1.append(temp_cytConc)
                if temp_conc > threshold_val:
                    conc1.append(temp_conc)      
                else:
                    conc1.append(np.nan)
           else:
                conc1.append(np.nan)
                cytConc1.append(np.nan)
           for voxels in range(int(center_voxel), int(center_voxel + width_voxel)):
               if an_y[t][voxels] > membrane_detect:
                  sum_mem2 = sum_mem2 + an_y[t][voxels]
                  sum_cyt2 = sum_cyt2 + Can_y[t][voxels]
                  spine2_voxels.append(voxels)
           if len(spine2_voxels) > 0:      
                temp_conc = (sum_cyt2 + sum_mem2) / (Na * v_voxel * len(spine2_voxels))
                temp_cytConc = (sum_cyt2) / (Na * v_voxel * len(spine2_voxels))
                cytConc2.append(temp_cytConc)
                if temp_conc > threshold_val:
                   conc2.append(temp_conc)      
                else:   
                   conc2.append(np.nan) 
           else:
                conc2.append(np.nan)
                cytConc2.append(np.nan)
           for x in range(int(center_voxel - width_voxel), int(center_voxel + width_voxel)):
               time_heat.append(round(chemTimes[t], 2))
               x_heat.append(round(x_axis[x], 2))
               if an_y[t][x] > membrane_detect:
                   #value_heat.append(1e3 * an_y[t][x] / (v_voxel * Na))
                   if x < center_voxel:
                       if sum_mem1 > 0:
                           value_heat.append(sum_mem1 / len(spine1_voxels))
                   else:
                       if sum_mem2 > 0:
                           value_heat.append(sum_mem2 / len(spine2_voxels))
               else:
                   value_heat.append(np.nan)
                   
       df["x"] = x_heat
       df["time"] = time_heat
       df["value"] = value_heat
       print(df)
       print(center_voxel)
       print(width_voxel)
       return df, recorded_chemTimes, conc1, conc2, cytConc1, cytConc2


def extract_shapes(cytConc, N, spine_create_delay, spine_spacing):
        tag_string = "CytConc_" + str(cytConc * 1e-3) + "_ns_" + str(N) + "_spine_create_delay_" + str(spine_create_delay) + "_spine_spacing_" + str(spine_spacing)
        phi = pd.read_csv("./" + tag_string + "phi.csv")
        hmean = pd.read_csv("./" + tag_string + "hmean.csv")
        theta = pd.read_csv("./" + tag_string + "theta.csv")
        saddle_start = pd.read_csv(tag_string + "sstart.csv")
        cConc = pd.read_csv(tag_string + "cConc.csv")
        rp = pd.read_csv(tag_string + "rp.csv")
        shapeTimes = pd.read_csv(tag_string + "shapeTimes.csv")
        x_save = []
        y_save = []
        times = []
        for i in range(len(phi)):
           saddle_start_list = remove_nan(saddle_start.iloc[i][1:])
           theta_list = remove_nan(theta.iloc[i][1:]) 
           Hmean_list = remove_nan(hmean.iloc[i][1:])
           Rp_list = remove_nan(rp.iloc[i][1:])
           x, y = dendShape.get_shape(saddle_start_list, theta_list, Hmean_list, Rp_list)
           non_zero = np.where(y > 1e-2)
           x_save.append(x[non_zero])
           y_save.append(y[non_zero])
           times.append(shapeTimes.iloc[i]["0"])
        return x_save, y_save, times


def heatmap(list1, list2, list1vr, list2vr, constvr, constvr_names, directory):
    df = pd.DataFrame()
    output_list = []
    vr1 = []
    vr2 = []
    num_list = []
    phitot_list = []
    energy_list = []
    ag_energy_list = []
    mismatch_energy_list = []
    entropy_energy_list = []
    for cvr in range(len(constvr_names)):
        globals()[constvr_names[cvr]] = constvr[cvr]
    for cvr in range(len(constvr_names)):
        print(constvr_names[cvr], globals()[constvr_names[cvr]])
    for l1 in range(len(list1)):    
        for l2 in range(len(list2)):    
             globals()[list1vr] = list1[l1]
             globals()[list2vr] = list2[l2]
             print(list1vr)
             print(list2vr)
             print(globals()[list1vr])
             print(globals()[list2vr])
             tag_string = "CytConc_" + str(round(Conc * 1e-3, 5)) + "_ns_" + str(N) + "_spine_create_delay_" + str(spine_create_delay) + "_spine_spacing_" + str(spine_spacing)

             vr1.append(globals()[list1vr])
             vr2.append(globals()[list2vr])
             print(tag_string) 
             df_shapeTime = pd.read_csv("./" + directory + tag_string + "shapeTimes.csv")
             try:
                 #plot_index = df_shapeTime[abs(df_shapeTime["0"] - time_of_plot) < 1].index[0]
                 plot_index = -1
             except:
                 print("Spine didn't form")
             phi = pd.read_csv("./" + directory + tag_string + "phi.csv")
             hmean = pd.read_csv("./" + directory + tag_string + "hmean.csv")
             theta = pd.read_csv("./" + directory + tag_string + "theta.csv")
             saddle_start = pd.read_csv("./" + directory + tag_string + "sstart.csv")
             cConc = pd.read_csv("./" + directory + tag_string + "cConc.csv")
             try:
               rp = pd.read_csv("./" + directory + tag_string + "rp.csv")
               last_phi = remove_nan(phi.iloc[plot_index][1:])
               last_hmean = remove_nan(hmean.iloc[plot_index][1:])
               last_rd = 1 / np.asarray(last_hmean)
               last_rp = remove_nan(rp.iloc[plot_index][1:])
               last_theta = remove_nan(theta.iloc[plot_index][1:])
               last_saddle_start = remove_nan(saddle_start.iloc[plot_index][1:])
               last_conc = remove_nan(cConc.iloc[plot_index][1:])
               spine_number = len(last_phi)
               num_list.append(spine_number)
               area = 2 * np.pi * last_rd**2 * (1 - np.cos(last_theta))
               print(area)
               phitot = np.asarray(last_phi) * phisat * np.asarray(area)
               phitot_list.append(sum(phitot))
               sum_N_energy = 0
               sum_ag_energy = 0
               sum_mismatch_energy = 0
               sum_entropy_energy = 0
               print("Num spines: ", spine_number)
               for n in range(spine_number):
                   sum_N_energy = sum_N_energy + EP.total_energy([ last_rd[n] * 1e6 * last_theta[n], last_theta[n], last_rp[n] * 1e6 ], phitot[n], last_conc[n], -45.0, -4.5)
                   sum_ag_energy = sum_ag_energy + EP.aggregation_energy(last_phi[n], phitot[n], area[n], -45.0) / EP.KBT
                   sum_mismatch_energy = sum_mismatch_energy + EP.mismatch(last_rd[n], last_phi[n], phitot[n], area[n]) / EP.KBT 
                   sum_entropy_energy = sum_entropy_energy + EP.entropy(last_phi[n], phitot[n], area[n], EP.khat) / EP.KBT

               energy_list.append(sum_N_energy)    
               ag_energy_list.append(sum_ag_energy)    
               mismatch_energy_list.append(sum_mismatch_energy)
               entropy_energy_list.append(sum_entropy_energy)
               
             except:  
               raise  
               print("No spines")
               energy_list.append(np.nan)
               print(plot_index)
               num_list.append(0)
               phitot_list.append(np.nan)
    df[list1vr] = vr1     
    df[list2vr] = vr2
    df["spines"] = num_list
    df["phitot"] = phitot_list
    df["Energy"] = energy_list
    df["Agg. Energy"] = ag_energy_list
    df["Entropy"] = entropy_energy_list
    df["Mismatch Energy"] = mismatch_energy_list
    df.to_csv("./data.csv")
    return df

figH = plt.figure(figsize = (18, 24))
ax1H = plt.subplot2grid((3, 3), (0, 0), colspan = 1) 
ax11H = plt.subplot2grid((3, 3), (0, 1), colspan = 1) 
ax12H = plt.subplot2grid((3, 3), (0, 2), colspan = 1) 
ax21H = plt.subplot2grid((3, 3), (1, 0), rowspan = 1) 
ax22H = plt.subplot2grid((3, 3), (1, 1), rowspan = 1) 
ax23H = plt.subplot2grid((3, 3), (1, 2), rowspan = 1) 
ax31H = plt.subplot2grid((3, 3), (2, 0), rowspan = 1) 
ax32H = plt.subplot2grid((3, 3), (2, 1), rowspan = 1) 
ax33H = plt.subplot2grid((3, 3), (2, 2), rowspan = 1) 
#ax41H = plt.subplot2grid((4, 3), (3, 0), rowspan = 2) 
#ax42H = plt.subplot2grid((4, 3), (3, 1), rowspan = 1) 
#ax43H = plt.subplot2grid((4, 3), (3, 2), rowspan = 1) 

ax1H.text(-0.2, 1.1, 'A', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax1H.transAxes)

ax11H.text(-0.2, 1.1, 'B', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax11H.transAxes)

ax12H.text(-0.2, 1.1, 'C', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax12H.transAxes)

ax21H.text(-0.2, 1.1, 'D', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax21H.transAxes)

ax22H.text(-0.2, 1.1, 'E', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax22H.transAxes)

ax23H.text(-0.25, 1.1, 'F', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax23H.transAxes)

ax31H.text(-0.2, 1.1, 'G', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax31H.transAxes)

ax32H.text(-0.2, 1.1, 'H', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax32H.transAxes)

ax33H.text(-0.2, 1.1, 'I',
        horizontalalignment='left', weight = "bold",
        verticalalignment='bottom',
        transform=ax33H.transAxes)


ax11H.set_xticks(cyt_conc_list)  
#def extract_conc_profile(cytConc, N, spine_create_delay, spine_spacing, filename):
def plot_conc_heatmap(conc, N, spine_create_delay, spine_spacing, molecule, axis, center, width, zoom_time, ticklabels, bar_vis, to_plot, cbar_flag):
    df_conc, chemTimes, conc1, conc2, cytConc1, cytConc2 = extract_conc_profile(conc, N, spine_create_delay, spine_spacing, molecule, center, width, zoom_time)       
    print("Lengths")
    print(len(df_conc["time"]))
    print(len(df_conc["x"]))
    print(len(df_conc["value"]))
    print(df_conc["x"])
    res = df_conc.pivot(index = "time", columns = "x", values = "value")
    print(res)
    if to_plot:
       #if cbar_flag == False:
       print("CBAR FLAG: ", cbar_flag)
       axis = sns.heatmap(res, ax = axis, cmap = 'cividis', cbar = cbar_flag, vmin = 0, vmax = 100, cbar_kws = {"label" : "#Membrane molecules"})
       #axis = sns.heatmap(res, ax = axis, cmap = 'cividis', cbar = bar_vis)
       if ticklabels:
          axis.set_xticks([0.5 * width * 1e-6 / diffL, int((width) * 1e-6 / diffL), (2 * width - 0.5 * width) * 1e-6 / diffL])
          axis.set_yticks([0, 0.5 * len(chemTimes), len(chemTimes)])
          axis.set_yticklabels([0, math.ceil(0.5 * max(df_conc["time"])), math.ceil(max(df_conc["time"]))])
          axis.set_xticklabels([center - 0.5 * width, center, center + 0.5 * width])
          axis.invert_yaxis()
          #axis.set_title(str(molecule) + " conc.)")
          axis.set_ylabel("t $(s)$")
          axis.set_xlabel("x ($\mu m$)")
       else:  
          axis.set_xticks([])  
          axis.set_yticks([])
          axis.set_ylabel("")
          axis.set_xlabel("")
    return chemTimes, conc1, conc2, cytConc1, cytConc2, df_conc   

def plot_shapes(conc, N, spine_create_delay, spine_spacing, axis):
   x_save, y_save, time = extract_shapes(conc, N, spine_create_delay, spine_spacing)
   maxH = max(y_save[-1])
   for i in range(len(x_save)):
      if len(y_save[i]) > 0:
         alpha = max(y_save[i]) / (maxH + 0.5)
         axis.plot(x_save[i], y_save[i], alpha = alpha, color = 'r')
   axis.set_xlabel("x")
   axis.set_ylabel("Height $\mu m$") 





def energy_compare_N():
   energy_compare_conc_list = [0.4]
   #def extract_data(mu0, k_agg1, scale_length, conc):
   N4eccl_energy_list = []
   N2eccl_energy_list = []
   N1eccl_energy_list = []
   select_eccl_shape = 0.4
   ax11H.set_ylim([0, 5])
   ax11H.set_xlim([4, 7])
   ax11H.set_xticks([])
   ax11H.set_yticks([])
   dend_plot_height = 0.5
   shape_lists_x = []
   shape_lists_y = []
   for eccl in energy_compare_conc_list:
       N4E, HRM, HTHETA, HRP, LRM, LTHETA, LRP = NCE.extract_data(-4.5, -60.0, 4, eccl)
       if eccl == select_eccl_shape:
         x, y, dms, dme = dendShape.dendShape( 0, HTHETA, HRM, HRP )
         axin11H_1 = ax11H.inset_axes([4, 4.0, 3.0, 1.0], transform = ax11H.transData) 
         axin11H_1.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
         #axin11H_1.add_patch( Rectangle( (0, 0), 10, dend_plot_height, fc = 'none', ec = 'grey') )
         plot_offset = 4.0
         shape_lists_x.append(x)
         shape_lists_y.append(y)
         for nplot in range(4):
            axin11H_1.plot(plot_offset + np.asarray(x) * 1e6, np.asarray(y) * 1e6 + dend_plot_height, "b", label = str(4), color = "b")
            plot_offset = plot_offset + 0.4
         #x, y, dms, dme = dendShape.dendShape( 0, LTHETA, LRM, LRP )
         axin11H_1.set_ylim([0, dend_plot_height + 0.15])
         axin11H_1.set_xlim([4, 6.0])
         axin11H_1.set_xticks([])
         axin11H_1.set_yticks([])
         #axin11H_1.plot([0, 10], [0] * 2, 'k', linestyle = "-.", linewidth = 6) 
         #ax22H.plot(np.asarray(x) * 1e6, np.asarray(y) * 1e6, "k")
       N4eccl_energy_list.append( N4E )
       N2E, HRM, HTHETA, HRP, LRM, LTHETA, LRP = NCE.extract_data(-4.5, -60.0, 2, eccl)
       if eccl == select_eccl_shape:
         x, y, dms, dme = dendShape.dendShape( 0, HTHETA, HRM, HRP )
         axin11H_2 = ax11H.inset_axes([4, 2.0, 3.0, 1.0], transform = ax11H.transData) 
         axin11H_2.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
         #axin11H_2.add_patch( Rectangle( (0, 0), 10, dend_plot_height, fc = 'none', ec = 'grey') )
         plot_offset = 4
         shape_lists_x.append(x)
         shape_lists_y.append(y)
         for nplot in range(2):
            axin11H_2.plot(plot_offset + np.asarray(x) * 1e6, np.asarray(y) * 1e6 + dend_plot_height, "b", label = str(4), color = "orange")
            plot_offset = plot_offset + 0.7
         #x, y, dms, dme = dendShape.dendShape( 0, LTHETA, LRM, LRP )
         axin11H_2.set_ylim([0, dend_plot_height + 0.15])
         axin11H_2.set_xlim([4, 5.2])
         axin11H_2.set_xticks([])
         axin11H_2.set_yticks([])
         #axin11H_2.plot([0, 10], [0] * 2, 'k', linestyle = "-.", linewidth = 6) 
         #ax22H.plot(np.asarray(x) * 1e6, np.asarray(y) * 1e6, "b")
       N2eccl_energy_list.append( N2E )
       N1E, HRM, HTHETA, HRP, LRM, LTHETA, LRP = NCE.extract_data(-4.5, -60.0, 1, eccl)
       if eccl == select_eccl_shape:
         x, y, dms, dme = dendShape.dendShape( 0, HTHETA, HRM, HRP )
         axin11H_3 = ax11H.inset_axes([4, 0.0, 3.0, 1.0], transform = ax11H.transData) 
         axin11H_3.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
         #axin11H_3.add_patch( Rectangle( (0, 0), 10, dend_plot_height, fc = 'none', ec = 'grey') )
         plot_offset = 4
         shape_lists_x.append(x)
         shape_lists_y.append(y)
         for nplot in range(1):
            axin11H_3.plot(plot_offset + np.asarray(x) * 1e6, np.asarray(y) * 1e6 + dend_plot_height, "b", label = str(4), color = "g")
            plot_offset = plot_offset + 0.7
         #ax21H.plot(np.asarray(x) * 1e6, np.asarray(y) * 1e6, "r", label = str(1))
         axin11H_3.set_ylim([0, dend_plot_height + 0.3])
         axin11H_3.set_xlim([4, 5.2])
         axin11H_3.set_xticks([])
         axin11H_3.set_yticks([])
         x_array = np.asarray([0, 0.2e-6]) * 1e6 + 4.8
         print("Xarray: ", x_array)
         print("max y:", max(y)*1e6)
         axin11H_3.plot(x_array, [y[0] * 1e6 + 0.6] * 2, 'k')  
         axin11H_3.plot([x_array[0]] * 2, [0.6, 0.8], 'k')  
         axin11H_3.text(x_array[-1], y[0] * 1e6 + 0.6, "0.2 $\mu m$")
         axin11H_3.text(x_array[0], y[0] * 1e6 + 0.8, "0.2 $\mu m$")
         #axin11H_3.plot([0, 10], [0] * 2, 'k', linestyle = "-.", linewidth = 6) 
         #x, y, dms, dme = dendShape.dendShape( 0, LTHETA, LRM, LRP )
         #ax22H.plot(np.asarray(x) * 1e6, np.asarray(y) * 1e6, "r")
       N1eccl_energy_list.append( N1E )
   num_clusters = [4, 2, 1]    
   for sl in range(len(shape_lists_x)):
       ax12H.plot(np.asarray(shape_lists_x[sl]) * 1e6, np.asarray(shape_lists_y[sl]) * 1e6, label = "N: " + str(num_clusters[0]))
       num_clusters.pop(0)
   ax12H.spines[['right', 'top']].set_visible(False)
   ax12H.set_ylabel("height $\mu m$")
   ax12H.set_xlabel("x $\mu m$")
   ax12H.legend(frameon = False, loc = "upper right")
       
   #ax21H.set_ylim([0, 0.1])
   #ax22H.set_ylim([0, 0.1])
   #ax21H.set_xlabel("$\mu m$")
   #ax21H.set_ylabel("$\mu m$")
   ax11H.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
   #ax22H.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
   #ax22H.set_xlabel("$\mu m$")
   #ax12H.set_title("Sim. Energy deviation %")
   #ax12H.set_xlabel("Conc $\mu M$")
   #ax12H.spines[['right', 'top']].set_visible(False)
   
# The following are for energy comparison for N spines from simulation
   #markers = ["X", "o", "v"]
   #ax1_colors = ["b", "m", "g"]
   #ax1H.plot(energy_compare_conc_list, np.asarray( N1eccl_energy_list ) / 1e3, ax1_colors[0],  marker = markers[0], markerfacecolor = "None", label = "$N_{calc}$ = 1")
   #ax1_colors.pop(0)
   #markers.pop(0)
   #ax1H.plot(energy_compare_conc_list, np.asarray( N2eccl_energy_list ) / 1e3, ax1_colors[0], marker = markers[0], markerfacecolor = "None", label = "$N_{calc}$ = 2")
   #ax1_colors.pop(0)
   #markers.pop(0)
   #ax1H.plot(energy_compare_conc_list, np.asarray( N4eccl_energy_list ) / 1e3, ax1_colors[0], marker = markers[0], markerfacecolor = "None", label = "$N_{calc}$ = 4")
   #ax1_colors.pop(0)
   #markers.pop(0)
   #ax12H.plot(energy_compare_conc_list, abs( (np.asarray( N1_sim ) - np.asarray(N1eccl_energy_list)) / np.asarray(N1eccl_energy_list)) * 100, "b", marker = "X", markerfacecolor = "None", label = "N = 1")
   #ax12H.plot(energy_compare_conc_list, abs( (np.asarray( N2_sim ) - np.asarray(N2eccl_energy_list)) / np.asarray(N2eccl_energy_list)) * 100, "m", marker = "o", markerfacecolor = "None", label = "N = 2")
   #ax12H.plot(energy_compare_conc_list, abs( (np.asarray(N4_sim) - np.asarray(N4eccl_energy_list)) / np.asarray(N4eccl_energy_list)) * 100, "g", marker = "v", markerfacecolor = "None", label = "N = 4")
   #ax12H.legend(frameon = False)   
   ax12H.legend(frameon = False)
   

df = heatmap([0.55, 0.6, 0.65, 0.7, 0.75, 1.0, 1.25, 1.5], [0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.6, 0.7, 0.8, 1.0, 1.4, 1.8, 2.0], "spine_spacing", "Conc", [4, 0.0], ["N", "spine_create_delay"], "./")
res = df.pivot(index = "spine_spacing", columns = "Conc", values = "spines")
print("N4: ", res)
#myColors = ["#000000", "#808000", "#ffff00","#0000ff", "#00FFFF"]
myColors = ["#808000", "#00FF00","#0000ff", "#FFFF00"]
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
ax1H = sns.heatmap(res, ax = ax1H, cmap = cmap, cbar_kws={'ticks': [0, 1, 2, 3, 4]})
ax1H.set_yticklabels(ax1H.get_yticklabels(), rotation = 0)
ax1H.invert_yaxis()
ax1H.set_title("Final N")


energy_compare_N()

def axis_mod(axis):
   #axis.spines[['right', 'top', 'left', 'bottom']].set_visible(False)
   axis.spines[['right', 'left', 'top', 'bottom']].set_visible(False)
   axis.set_yticks([])
   axis.set_xticks([])

def non_zero_shape(x, y):
   x_nz = []
   y_nz = []
   for i in range(len(y)):
      if y[i] > 0:
        y_nz.append(y[i])
        x_nz.append(x[i])
   return x_nz, y_nz

def inset_shapes(axis, x0, y0, x1, y1, ref_x, ref_y, color):
   heatmap_inset = axis.inset_axes([ref_x, ref_y, 0.25, 0.25])
   heatmap_inset_twin = heatmap_inset.twinx()
   x0_nz, y0_nz = non_zero_shape(x0, y0)
   heatmap_inset.plot(x0_nz, 10 + np.asarray(y0_nz) * 1e3, color) 
   maxy0 = max(y0) * 1e3
   heatmap_inset.set_yticks([0, np.floor(maxy0)]) 
   heatmap_inset.set_xlim(8.5, 11.5)
   """
   if len(x1) > 0:
      heatmap_inset_twin.plot(x1, np.asarray(y1) * 1e3, color) 
      maxy1 = max(y1) * 1e3
      heatmap_inset_twin.set_yticks([0, maxy1])
   """
   x1_nz, y1_nz = non_zero_shape(x1, y1)
   heatmap_inset.plot(x1_nz, 10 + np.asarray(y1_nz) * 1e3, color) 
   axis_mod(heatmap_inset)
   axis_mod(heatmap_inset_twin)
   return heatmap_inset
   


def ref_line(time_ref, color, x):
      mark_line_offset = 8.5
      heatmap_axis.plot([ (x - mark_line_offset) / (diffL * 1e6), (x + 0.5 - mark_line_offset) / (diffL * 1e6)], [time_ref / 0.05 , time_ref/ 0.05 ], color = color, linestyle = "--")


def contiguous_chunk(df):
    x_1 = []
    x_2 = []
    first_half = df[df["x"] < 10.0]
    second_half = df[df["x"] > 10.0]
    print("FIRST HALF", first_half)
    print("SECOND HALF", second_half)
    for i in range(len(first_half)):
        if np.isnan(first_half.iloc[i]["value"]) == False:
           x_1.append(first_half.iloc[i]["x"])
    for i in range(len(second_half)):
        if np.isnan(second_half.iloc[i]["value"]) == False:
           x_2.append(second_half.iloc[i]["x"])
    return x_1, x_2

def mem_plot(df, time, axis):
      df = df[ df["time"] >= time ]
      time_available = round(df.iloc[0]["time"], 2)
      df = df[df["time"] == time_available]
      xm1, xm2 = contiguous_chunk(df)
      axis.plot([xm1[0], xm1[-1]], [0, 0], "y", linewidth = 10)
      axis.plot([xm2[0], xm2[-1]], [0, 0], "y", linewidth = 10)


spacing_competition = [0.6, 0.7, 0.8, 1.0, 2.0]
conc_list = [0.6, 0.8, 1.4]
conc_competition = 1.4
delay = 25.0
ticklabels = True
select_spacings = [0.6, 1.0, 2.0]
select_axes = [ax21H, ax22H, ax23H]
heatmap_axis = ax21H
colorbar_count = 0
ax32H_plot_count = 0
for cc in conc_list:
  comp_sp_list = []
  cyt_init_2 = []
  for sp in spacing_competition:
    plot_heatmap = False
    print("COLORBAR: ", colorbar_count)
    if colorbar_count == 0:
       cbar_flag = True
    else:
       cbar_flag = False
    print("Current conc and sp: ", cc, sp)
    if sp in select_spacings and cc == conc_competition:
       print("Plotting heatmap for membrane")
       plot_heatmap = True
       heatmap_axis = select_axes[0]
       select_axes.pop(0)
       select_spacings.pop(0)
       heatmap_axis.text(0.1, 0.99, str(cc) + ": $\mu M$," + str(sp) + ": $\mu m$", weight = "normal",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=heatmap_axis.transAxes)
     
       #heatmap_axis.set_title(str(cc) + ": $\mu M$," + str(sp) + ": $\mu m$")
       heatmap_axis.set_ylabel("t $(s)$")
       heatmap_axis.set_xlabel("x $(\mu m)$")
       print("Select axis: ", heatmap_axis)
       colorbar_count = colorbar_count + 1
    comp_sp_list.append(sp)
    if cc == conc_competition and sp == 2.0:
      x0, y0, time_2_init = get_parameters(cc * 1e-3, sp, 0, delay)
      x1, y1, time_2_init = get_parameters(cc * 1e-3, sp, 1, delay)
      heatmap_inset1 = inset_shapes(heatmap_axis, x0, y0, x1, y1, 1.0, 0.45, "b")
      ref_line(time_2_init, "b", 11.0)


      x0, y0 = get_parameters_select_time(cc * 1e-3, sp, time_2_init + 25, 0, delay)
      x1, y1 = get_parameters_select_time(cc * 1e-3, sp, time_2_init + 25, 1, delay)
      heatmap_inset2 = inset_shapes(heatmap_axis, x0, y0, x1, y1, 1.0, 0.85, "m")
      ref_line(time_2_init + 25, "m", 11.2)


      chemTimes, conc1, conc2, cytConc1, cytConc2, df_conc = plot_conc_heatmap(cc, 2, delay, sp, "membrane", heatmap_axis, 10, 1.5, [0, 60], ticklabels, True, plot_heatmap, cbar_flag)
      mem_plot(df_conc, time_2_init, heatmap_inset1)
      mem_plot(df_conc, time_2_init + 25, heatmap_inset2)

      axis_mod(heatmap_inset1)
    else:
      chemTimes, conc1, conc2, cytConc1, cytConc2, df_conc = plot_conc_heatmap(cc, 2, delay, sp, "membrane", heatmap_axis, 10, 1.5, [0, 60], ticklabels, True, plot_heatmap, cbar_flag)
    nonZ = np.where(np.asarray(conc2) > 0)[0]
    cyt_init_2.append(conc2[nonZ[0]] * 1e3)
  if ax32H_plot_count == 0:
     ax32H.plot(comp_sp_list, np.asarray(cyt_init_2), label = "Total Conc. \n" + str(cc) + "$\mu M$")
  else:
     ax32H.plot(comp_sp_list, np.asarray(cyt_init_2), label = str(cc) + "$\mu M$")
  ax32H_plot_count += 1
  
ax22H.set_ylabel("")
ax22H.set_yticks([])
ax23H.set_ylabel("")
ax23H.set_yticks([])
ax32H.legend()
ax32H.set_ylim(0, 0.8)
ax32H.set_ylabel("Conc. at spine2 $\mu M$")
ax32H.set_xlabel("Distance from spine1 $\mu m$")
ax32H.spines[['right', 'top']].set_visible(False)
ax32H.legend(frameon = False)
cc = 1.4
sp = 0.7
delay = 0.0
heatmap_axis = ax31H
plot_heatmap = True
x0, y0, time_2_merge = get_parameters_merge(cc * 1e-3, sp, 0, delay)
x1, y1, time_2_merge = get_parameters_merge(cc * 1e-3, sp, 1, delay)
heatmap_inset1 = inset_shapes(heatmap_axis, x0, y0, x1, y1, 0.85, 0.35, "b")
ref_line(time_2_merge, "b", 10.6)

x0, y0 = get_parameters_select_time(cc * 1e-3, sp, time_2_merge + 25, 0, delay)
#x1, y1 = get_parameters_select_time(cc * 1e-3, sp, time_2_merge + 10, 1)
heatmap_inset2 = inset_shapes(heatmap_axis, x0, y0, [], [], 0.85, 0.7, "m")
ref_line(time_2_merge + 20, "m", 10.6)

chemTimes, conc1, conc2, cytConc1, cytConc2, df_conc = plot_conc_heatmap(cc, 2, delay, sp, "membrane", heatmap_axis, 10, 1.5, [0, 60], ticklabels, False, plot_heatmap, cbar_flag)
mem_plot(df_conc, time_2_merge, heatmap_inset1)
mem_plot(df_conc, time_2_merge + 20, heatmap_inset2)
ax31H.text(0.0, 0.99, str(cc) + ": $\mu M$," + str(sp) + ": $\mu m$, \n delay = " + str(delay) + " s", weight = "normal",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax31H.transAxes)
#heatmap_axis.set_title(str(cc) + ": $\mu M$," + str(sp) + ": $\mu m$")
heatmap_axis.set_ylabel("t $(s)$")
heatmap_axis.set_xlabel("x $(\mu m)$")
print("CytInit: ", cyt_init_2)
def plot_exp_data():
   exp_data = pd.read_csv("./dist_data_LTP.csv")
   ax_histx = ax32H.inset_axes([0, 1.05, 1, 0.25], sharex=ax32H)
   ax_histy = ax32H.inset_axes([1.05, 0, 0.25, 1], sharey=ax32H)
   scatter_hist(exp_data["ISD"], exp_data["SD"], ax32H, ax_histx, ax_histy)
   ax32H.set_ylim([0, 2])
   ax32H.set_xlim([0, 2])
   ax32H.set_xlabel("ISD")
   ax32H.set_ylabel("SD")
   ax32H.set_xticks([0, 0.5, 1.0, 1.5])
   ax32H.set_yticks([0, 0.5, 1.0, 1.5])


def plot_axis(x_length, y_length):
    axis = ax2.inset_axes([-0.15, 0, 0.15, 0.15])
    axis.plot([0, x_length], [0, 0 * 1e-3], 'k')
    axis.plot([0, 0], [0, y_length], 'k')
    axis.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xlim([0, 20])
    axis.set_ylim([0, 55])
    #axis.set_xlabel(str(x_length) + "$\mu m$")
    axis.set_ylabel(str(y_length) + "nm")

def list_not_list(list1, list2): 
    print(list1)
    print(list2)
    for l1 in list1:
        if l1 not in list2:
           deleted_one = l1
           print("Deleted: ", deleted_one)
           break
    return deleted_one   

def deletion_times(tag_string):
    phi = pd.read_csv("./" + tag_string + "phi.csv")
    shapeTimes = pd.read_csv("./" + tag_string + "shapeTimes.csv")
    locations = pd.read_csv("./" + tag_string + "locations.csv")
    ordered_locations = np.asarray(locations.iloc[0][1:]) * 1e6
    ordered_location_indices = list(np.arange(0, num_spines, 1))
    max_spines = num_spines
    spines_present = []
    deletion_data = []
    for n in range(len(phi)):
        len_phi = len(remove_nan(phi.iloc[n][1:]))
        if len_phi < max_spines:
            new_locations = np.asarray(locations.iloc[n][1:]) * 1e6
            new_locations_indices = []
            for nl in new_locations:
                ordered_count = 0
                for ol in ordered_locations:
                    if nl == ol:
                        new_locations_indices.append(ordered_count)
                    ordered_count = ordered_count + 1    
            max_spines = len_phi
            try:
                deleted_one = list_not_list(temp, new_locations_indices)
                print("Deleted: ", deleted_one)
                deletion_data.append([deleted_one, shapeTimes.iloc[n]["0"]])
            except:
                print("Nothing")
            temp = new_locations_indices
            spines_present.append(new_locations_indices)
    return deletion_data
ax21H.spines[['left', 'bottom']].set_visible(True)
ax22H.spines[['left', 'bottom']].set_visible(True)
ax23H.spines[['left', 'bottom']].set_visible(True)
ax31H.spines[['left', 'bottom']].set_visible(True)
plt.subplots_adjust(hspace = 0.5, wspace = 0.5)



df_spacing = min_space(ax33H)

#% start: automatic generated code from pylustrator
plt.figure(1).ax_dict = {ax.get_label(): ax for ax in plt.figure(1).axes}
import matplotlib as mpl
getattr(plt.figure(1), '_pylustrator_init', lambda: ...)()
plt.figure(1).ax_dict["<colorbar>"].set(position=[0.8709, 0.07571, 0.01291, 0.1925])
plt.figure(1).axes[2].legend(loc=(0.6501, 0.8595), frameon=False)
plt.figure(1).axes[6].set(position=[0.125, 0.07681, 0.1937, 0.1925])
plt.figure(1).axes[7].set(position=[0.4156, 0.07681, 0.1937, 0.1925])
plt.figure(1).axes[8].set(position=[0.7063, 0.07681, 0.155, 0.1925])
#% end: automatic generated code from pylustrator
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import readXML as RX
import xml.etree.ElementTree as ET
import seaborn as sns
import energy_process as EP
import shankReac as sReac
import plot_shape
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
typeStim = "b"
time_of_plot = 60


def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

def heatmap(list1, list2, list1vr, list2vr, constvr, constvr_names, directory):
    df = pd.DataFrame()
    output_list = []
    vr1 = []
    vr2 = []
    num_list = []
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
             nb_per_loc =  int(globals()["nb"] / globals()["NstimL"])
             print("NB: ", nb_per_loc)
             spacing = globals()["spacing"]
             spine_spacing = globals()["spine_spacing"]
             if coactive == 0.0:
                 spine_spacing = spacing
             tag_string = "NstimL_" + str(NstimL + 1) + "_nb_" + str(nb) + "_spine_spacing_" + str(spine_spacing) + "_spacing_" + str(spacing) + "_freq_" + str(freq) + "_ypos_" + str(ypos)+ "_coactive_" + str(coactive)
             vr1.append(globals()[list1vr])
             vr2.append(globals()[list2vr])
             print(tag_string) 
             df_shapeTime = pd.read_csv("./" + directory + "/" + tag_string + "shapeTimes.csv")
             try:
                 #plot_index = df_shapeTime[abs(df_shapeTime["0"] - time_of_plot) < 1].index[0]
                 plot_index = -1
             except:
                 print("Spine didn't form")
             phi = pd.read_csv("./" + directory + "/" + tag_string + "phi.csv")
             hmean = pd.read_csv("./" + directory + "/" +  tag_string + "hmean.csv")
             theta = pd.read_csv("./" + directory + "/" + tag_string + "theta.csv")
             saddle_start = pd.read_csv("./" + directory + "/" + tag_string + "sstart.csv")
             try:
               rp = pd.read_csv("./" + directory + "/" + tag_string + "rp.csv")
               last_phi = remove_nan(phi.iloc[plot_index][1:])
               last_hmean = remove_nan(hmean.iloc[plot_index][1:])
               last_rd = 1 / np.asarray(last_hmean)
               last_rp = remove_nan(rp.iloc[plot_index][1:])
               last_theta = remove_nan(theta.iloc[plot_index][1:])
               last_saddle_start = remove_nan(saddle_start.iloc[plot_index][1:])
               spine_number = len(last_phi)
               num_list.append(spine_number)
             except:  
               print("No spines")
               print(plot_index)
               num_list.append(0)
    df[list1vr] = vr1     
    df[list2vr] = vr2
    df["spines"] = num_list
    print(df)
    return df

figH = plt.figure(figsize = (10, 8))
ax1H = plt.subplot2grid((2, 3), (0, 0), colspan = 1) 
ax11H = plt.subplot2grid((2, 3), (0, 1), colspan = 1) 
ax12H = plt.subplot2grid((2, 3), (0, 2), colspan = 1) 
ax2H = plt.subplot2grid((2, 3), (1, 0), colspan = 2) 

ax1H.text(-0.2, 1.1, 'A',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax1H.transAxes)



df = heatmap([0.4, 0.5, 0.6, 0.75, 1.0, 2.0], [0.75, 1.0, 2.0, 4.0, 6.0, 8.0], "freq", "spacing", [8, 2, 1.0, 3, 0.25, 0.0, 0.75], ["nb", "num_spines", "NstimL", "TK", "ypos", "coactive", "spine_spacing"], "./")
res = df.pivot(index = "freq", columns = "spacing", values = "spines")
ax1H = sns.heatmap(res, ax = ax1H, cmap = 'cividis', annot = True, cbar_kws={'ticks': [1, 2]}, xticklabels = [0.75, 1.0, 2.0, 4.0, 6.0, 8.0], yticklabels = [0.4, 0.5, 0.6, 0.75, 1.0, 2.0])
ax1H.set_xticklabels(ax1H.get_xticklabels(), rotation = 90)
ax1H.set_yticklabels(ax1H.get_yticklabels(), rotation = 0)
ax1H.invert_yaxis()

print("Plotting  coactive")
df = heatmap([2.0, 3.0, 4.0], [0.3, 0.5, 0.7, 1.0, 2.0, 2.5, 3.0], "NstimL", "spacing", [1, 2, 3, 0.25, 1.0, 0.75, 0.25], ["nb", "num_spines", "TK", "ypos", "coactive", "spine_spacing", "freq"], "./")
res = df.pivot(index = "NstimL", columns = "spacing", values = "spines")
ax11H = sns.heatmap(res, ax = ax11H, cmap = 'cividis', annot = True, cbar_kws={'ticks': [1, 2]}, yticklabels = [2.0, 3.0, 4.0], xticklabels = [0.3, 0.5, 0.7, 1.0, 2.0, 2.5, 3.0])
ax11H.set_xticklabels(ax11H.get_xticklabels(), rotation = 90)
ax11H.set_yticklabels(ax11H.get_yticklabels(), rotation = 0)
ax11H.invert_yaxis()
"""
df = heatmap([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], [1.0], "NstimL", "spacing", [1, 15, 3, 0.25, 1.0, 0.75, 0.25], ["nb", "num_spines", "TK", "ypos", "coactive", "spine_spacing", "freq"], "../chemistry_min_spacing_maxSpines")
ax12H.plot(df["NstimL"], df["spines"])
ax12H.set_xlabel("NstimL")
ax12H.set_ylabel("N")
ax12H.spines[['right', 'top']].set_visible(False)
ax12H.set_xticks([2, 3, 4, 5, 6, 7, 8, 9])
"""
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

plt.show()

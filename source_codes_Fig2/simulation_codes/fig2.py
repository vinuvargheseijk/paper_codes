import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import totalEnergy as TE
import numpy as np
import bisect
import writeXML
import readXML as RX
import xml.etree.ElementTree as ET
#import pylustrator
#pylustrator.start()

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import get_min_phitot as GM 
sns.set_context('poster')
select_kagg = -60.0
figH = plt.figure(figsize = (11.4, 24))
ax11 = plt.subplot2grid((4, 2), (0, 0), colspan = 1)
ax12 = plt.subplot2grid((4, 2), (0, 1), colspan = 1)
ax21 = plt.subplot2grid((4, 2), (1, 0), colspan = 1)
ax22 = plt.subplot2grid((4, 2), (1, 1), rowspan = 1)
ax31 = plt.subplot2grid((4, 2), (2, 0), rowspan = 1)
ax32 = plt.subplot2grid((4, 2), (2, 1), rowspan = 1)
ax41 = plt.subplot2grid((4, 2), (3, 0), rowspan = 1)
ax42 = plt.subplot2grid((4, 2), (3, 1), rowspan = 1)

ax11.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax12.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax21.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax22.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax11.set_xticks([])
ax12.set_xticks([])
ax11.set_yticks([])
ax12.set_yticks([])

ax11.text(-0.1, 1.1, 'A', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax11.transAxes)

ax12.text(-0.1, 1.1, 'B', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax12.transAxes)

ax21.text(-0.1, 1.1, 'C', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax21.transAxes)

ax22.text(-0.1, 1.1, 'D', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax22.transAxes)

ax31.text(-0.1, 1.1, 'E', weight = "bold", 
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax31.transAxes)

ax32.text(-0.1, 1.1, 'F', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax32.transAxes)

ax41.text(-0.1, 1.1, 'G', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax41.transAxes)

ax42.text(-0.1, 1.1, 'H', weight = "bold", 
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax42.transAxes)

def compare_phitot(res_calc, res_sim):
    return round(abs(res_calc - res_sim), 1)

def line_plots(v_parameters, c_parameters, df, plot_name):
    plot_list = []
    for vp1 in v_parameters[1]:
        f = df[df[v_parameters[0]] == vp1]
        print(f)
        f = f[f[c_parameters[0]] == c_parameters[1]]
        print(f)
        try:
          print("PHITOT: ", f[plot_name].iloc[0])
          plot_list.append(f[plot_name].iloc[0])
        except:  
          plot_list.append(np.nan)
    return plot_list           


def line_plots_sim(v_parameters, c_parameters, phi_entire, df, plot_name):
    plot_list = []
    f = df[df["phi_entire"] == phi_entire]
    print(f)
    for vp1 in v_parameters[1]:
        g = f[f[v_parameters[0]] == vp1]
        print(g)
        h = g[g[c_parameters[0]] == c_parameters[1]]
        print(h)
        try:
           plot_list.append(h[plot_name].iloc[0])
        except:
           plot_list.append(np.nan)
    return plot_list

def pivot_df(df, index_v, columns_v, values_v):
    return df.pivot(index = index_v, columns = columns_v, values = values_v)

energy_calc_dir = "../"

def optimalPhitot_calc(phi_entire, kagg): #Extracting from energy landscape calculations
    file = "mu-4.5_probability_" + str(kagg) + "_" + str(phi_entire) + "_" + "5.0nm.csv"
    df = pd.read_csv(energy_calc_dir + file)
    energy_list = list(df["highEnergy"])
    minEnergy_index = energy_list.index(np.nanmin(energy_list))
    minEnergy = np.nanmin(energy_list)
    phitot_list = list(df["phitot"])
    optimalPhitot = phitot_list[minEnergy_index]
    return optimalPhitot, minEnergy



colors = ['k', 'r', 'b', 'g', 'c', 'm', 'y']
phi_entire_list = [2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0]

landscape_min_phitot = []
landscape_min_energy = []
for pe in phi_entire_list:
   o_phitot, o_energy = optimalPhitot_calc(pe, select_kagg)
   landscape_min_phitot.append(o_phitot)
   landscape_min_energy.append(o_energy)
   print("Optimal phitot and energy for: " + str(pe) +": ", o_phitot, o_energy)

agg_list = [-45.0, -50.0, -55.0, -60.0, -65.0]
mu0_list = [-4.5]
var_param1 = agg_list
var_param2 = phi_entire_list
const_param1 = mu0_list[0]
const_param_name = "mu0"

markers = ["o", "X", "*"]
error_plot = 1
ccount = 0
df_sim = pd.read_csv("./sim.csv")
wave_speed_demo_k_aggs = [-50.0, -60.0]
wave_speed_demo_phi_entire = [2000.0, 3000.0, 4000.0, 6000.0, 8000.0]
wave_speed_conc_list = 1e3 * np.asarray(wave_speed_demo_phi_entire) / (TE.V_cyt *  TE.Na)
dfWs1 = pd.read_csv("../sample_plot_wave_speed/" + str(wave_speed_demo_phi_entire[0]) + "_" + str(wave_speed_demo_k_aggs[1]) + ".csv")
dfWs2 = pd.read_csv("../sample_plot_wave_speed/" + str(wave_speed_demo_phi_entire[1]) + "_" + str(wave_speed_demo_k_aggs[1]) + ".csv")
dfWs3 = pd.read_csv("../sample_plot_wave_speed/" + str(wave_speed_demo_phi_entire[2]) + "_" + str(wave_speed_demo_k_aggs[1]) + ".csv")
dfWs4 = pd.read_csv("../sample_plot_wave_speed/" + str(wave_speed_demo_phi_entire[3]) + "_" + str(wave_speed_demo_k_aggs[1]) + ".csv")

axins211 = inset_axes(ax11, width="50%", height="75%", bbox_to_anchor=(0, .6, 1.8, .5), bbox_transform=ax11.transAxes, loc=3)
axins212 = inset_axes(ax11, width="50%", height="75%", bbox_to_anchor=(0, 0.0, 1.8, .5), bbox_transform=ax11.transAxes, loc=3)
ax11.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax11.set_yticks([])
ax11.set_xticks([])
ax11.set_xlabel("x $\mu m$", labelpad = 15)
ax11.set_ylabel("Conc $\mu M$", labelpad = 15)

Turing = "../abstract_plots_fig2/Turing_50_6e-05_1_1_5000.0.csv"
wave = "../abstract_plots_fig2/wave_50_6e-05_1_1.csv"

df_turing = pd.read_csv(Turing)
df_wave = pd.read_csv(wave)

wsIntensity = [1, 2, 3]
def arrow_plot(axis, mag_index, num_stim):
  axis.arrow(60, 10, 0, -2,
          head_width = 2 * wsIntensity[mag_index],
          width = 0.25 * wsIntensity[mag_index],
          fc = 'red',   
          ec ='red')
  if num_stim > 1:
      axis.arrow(80, 10, 0, -2,
          head_width = 2 * wsIntensity[0],
          width = 0.25 * wsIntensity[0],
          fc = 'red',   
          ec ='red')
      axis.arrow(100, 10, 0, -2,
          head_width = 2 * wsIntensity[0],
          width = 0.25 * wsIntensity[0],
          fc = 'red',   
          ec ='red')

def marker_stim(axis, mag_index, num_stim):
    axis.plot(60, 13, marker = "v", markerfacecolor = "red", markeredgecolor = "red", markersize = 8 * wsIntensity[mag_index])
    if num_stim > 1:
       axis.plot(80, 13, marker = "v", markerfacecolor = "red", markeredgecolor = "red", markersize = 8 * wsIntensity[0])
       axis.plot(100, 13, marker = "v", markerfacecolor = "red", markeredgecolor = "red", markersize = 8 * wsIntensity[0])
          


axins211.set_title("Turing")
axins211.plot(np.asarray(df_turing["x"]) * 1e6, df_turing["A"])
axins212.plot(np.asarray(df_wave["x"]) * 1e6, df_wave["A"], color = "k")
axins212.set_title("Wave pinning")
axins211.spines[['top', 'right']].set_visible(False)
axins212.spines[['top', 'right']].set_visible(False)
axins211.set_xticks([])

marker_stim(axins211, 0, 1)  
marker_stim(axins212, 0, 1)  

wave3 = "../abstract_plots_fig2/wave_50_6e-05_3_" + str(wsIntensity[0]) + ".csv"
wave32 = "../abstract_plots_fig2/wave_50_6e-05_3_" + str(wsIntensity[1]) + ".csv"
wave33 = "../abstract_plots_fig2/wave_50_6e-05_3_" + str(wsIntensity[2]) + ".csv"

df_wave3 = pd.read_csv(wave3)
df_wave32 = pd.read_csv(wave32)
df_wave33 = pd.read_csv(wave33)
axins121 = inset_axes(ax12, width="50%", height="75%", bbox_to_anchor=(0, 0.7, 1.8, .3), bbox_transform=ax12.transAxes, loc=3)
axins122 = inset_axes(ax12, width="50%", height="75%", bbox_to_anchor=(0, 0.35, 1.8, .3), bbox_transform=ax12.transAxes, loc=3)
axins123 = inset_axes(ax12, width="50%", height="75%", bbox_to_anchor=(0, 0.0, 1.8, .3), bbox_transform=ax12.transAxes, loc=3)
ax12.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
axins121.spines[['top', 'right']].set_visible(False)
axins122.spines[['top', 'right']].set_visible(False)
axins123.spines[['top', 'right']].set_visible(False)
axins121.set_title("Wave pinning")
ax12.set_yticks([])
ax12.set_xticks([])
axins121.plot(np.asarray(df_wave3["x"]) * 1e6, df_wave3["A"], 'k')
axins122.plot(np.asarray(df_wave32["x"]) * 1e6, df_wave32["A"], 'k')
axins123.plot(np.asarray(df_wave33["x"]) * 1e6, df_wave33["A"], 'k')
axins122.set_xticks([])
axins121.set_xticks([])
axins121.set_ylim([0.01, 15])
axins122.set_ylim([0.01, 15])
axins123.set_ylim([0.01, 15])

ax21.set_yticks([])
ax21.set_xticks([])
ax22.set_yticks([])
ax22.set_xticks([])
ax12.set_xlabel("x $\mu m$", labelpad = 15)


#arrow_plot(axins121, 0, 3)  
#arrow_plot(axins122, 1, 3)  
#arrow_plot(axins123, 2, 3)  
marker_stim(axins121, 0, 3)
marker_stim(axins122, 1, 3)
marker_stim(axins123, 2, 3)
  

def find_intersect(values, phitots):
    lst_values = list(values)
    lst_values.sort()
    index = bisect.bisect_left(lst_values, 0)
    phitot_index = np.where(values == lst_values[index])
    print(phitot_index)
    return index, values[index], phitots[phitot_index[0]]

directory = "../sample_shapes_merging"
filename_X = "CytConc_0.0014_ns_2_spine_create_delay_25.0_spine_spacing_0.8X.xml"
filename_Y = "CytConc_0.0014_ns_2_spine_create_delay_25.0_spine_spacing_0.8Y.xml"
file_saddle = "CytConc_0.0014_ns_2_spine_create_delay_25.0_spine_spacing_0.8sstart.csv"
df_saddle = pd.read_csv(directory + "/" + file_saddle)
treeX = ET.parse(directory + "/" + filename_X)
tree = [treeX]
sum_an,max_an_y,an_y=RX.plotXML(filename_X,tree)
X_anY = an_y
treeY = ET.parse(directory + "/" + filename_Y)
tree = [treeY]
sum_an,max_an_y,an_y=RX.plotXML(filename_Y,tree)
Y_anY = an_y
axins221 = inset_axes(ax22, width="50%", height="75%", bbox_to_anchor=(0, .6, 1.8, .5), bbox_transform=ax22.transAxes, loc=3)
axins222 = inset_axes(ax22, width="50%", height="75%", bbox_to_anchor=(0, 0.0, 1.8, .5), bbox_transform=ax22.transAxes, loc=3)
#axins223 = inset_axes(ax22, width="25%", height="25%", bbox_to_anchor=(0.6, 0.3, 1.8, .5), bbox_transform=ax22.transAxes, loc=3)

axins221.plot(X_anY[249], Y_anY[249])
axins222.plot(X_anY[-1], Y_anY[-1])
#axins223.plot(X_anY[249], Y_anY[249])
nonzeros = np.where(np.asarray(Y_anY[249]) > 0.001)
axins221.set_xlim([9, 11])
axins222.set_xlim([9, 11])
axins221.set_ylim([0.0005, 0.2])
axins222.set_ylim([0.0005, 0.2])
#axins223.set_ylim([0, 0.002])
axins221.plot([df_saddle.iloc[249]["1"] * 1e6] * 10, np.linspace(0, max(Y_anY[nonzeros[0][0]]), 10), 'k', linestyle = "-.")
axins221.spines[['top', 'right']].set_visible(False)
axins222.spines[['top', 'right']].set_visible(False)
#axins223.spines[['top', 'right', 'left']].set_visible(False)
axins221.set_xticks([])
axins221.set_yticks([0, 0.1])
axins222.set_yticks([0, 0.1])
axins222.set_xlabel("x $\mu m$")
axins222.set_ylabel("$\mu m$")
axins221.set_ylabel("$\mu m$")


ax31.plot(dfWs1["phitot"], np.asarray(dfWs1["Speed"]) * 1e9, label = str(round(wave_speed_conc_list[0],1)) + "$\mu M$")
intersect, intersect_speed, intersect_phitot = find_intersect(dfWs1["Speed"], dfWs1["phitot"])
ax31.plot(intersect_phitot, intersect_speed, marker = "o", markerfacecolor = "k", markeredgecolor = "k")
ax31.plot(dfWs2["phitot"], np.asarray(dfWs2["Speed"]) * 1e9, label = str(round(wave_speed_conc_list[1],1)))
intersect, intersect_speed, intersect_phitot = find_intersect(dfWs2["Speed"], dfWs2["phitot"])
ax31.plot(intersect_phitot, intersect_speed, marker = "o", markerfacecolor = "k", markeredgecolor = "k")
ax31.plot(dfWs3["phitot"], np.asarray(dfWs3["Speed"]) * 1e9, label = str(round(wave_speed_conc_list[2],1)))
intersect, intersect_speed, intersect_phitot = find_intersect(dfWs3["Speed"], dfWs3["phitot"])
ax31.plot(intersect_phitot, intersect_speed, marker = "o", markerfacecolor = "k", markeredgecolor = "k")
ax31.plot(dfWs4["phitot"], np.asarray(dfWs4["Speed"]) * 1e9, label = str(round(wave_speed_conc_list[3],1)))
intersect, intersect_speed, intersect_phitot = find_intersect(dfWs4["Speed"], dfWs4["phitot"])
ax31.plot(intersect_phitot, intersect_speed, marker = "o", markerfacecolor = "k", markeredgecolor = "k")

ax31.set_ylabel("$W_{s}$ nm/s")
ax31.set_xlabel("$\phi_{tot}$")
ax31.plot(dfWs4["phitot"], [0] * len(dfWs4["phitot"]), 'k', linestyle = "--")
ax31.spines[['top', 'right']].set_visible(False)
ax32.spines[['top', 'right']].set_visible(False)
ax31.legend(frameon = False)

line_plot = False
heatmap = True

#The following commented lines are fore plottinh heatmap from wave speed calcs
"""
for vp2 in var_param2:
  df = pd.read_csv("./Ws_" + str(vp2) + ".csv")
  print(df_sim)
  #df = df.drop_duplicates()
  cyt_conc = (vp2 / (TE.V_cyt * TE.Na)) * 1e3  
  color = colors[ccount]
  #line_plots(v_parameters, c_parameters, df, plot_name):
  plot_list = line_plots(["Agg. coeff.", var_param1], [const_param_name, const_param1], df, "Phitot")
  plot_list_Ws = line_plots(["Agg. coeff.", var_param1], [const_param_name, const_param1], df, "Speed")
  plot_list_size = line_plots_sim(["Aggregation", var_param1], [const_param_name, const_param1], vp2, df_sim, "protSize")
  print(plot_list_size)
  plot_list_sim = line_plots_sim(["Aggregation", var_param1], [const_param_name, const_param1], vp2, df_sim, "phitot")
  print("PLOT LIST: ", plot_list)
  print(cyt_conc)
  print("PLOT SIM LIST: ", plot_list_sim)
  #ax22.plot(agg_list, np.asarray(plot_list_Ws), marker = "o", label = "Conc. " + str(round(cyt_conc,1)))
  #ax22.spines[['top', 'right']].set_visible(False)

if heatmap == True:
    calc_conc_list = []
    calc_agg_list = []
    calc_speed_list = []
    calc_phitot_list = []
    calc_phi_entire_list = []
    calc_Energy_list = []
    for pe in phi_entire_list:
      df = pd.read_csv("./Ws_" + str(pe) + ".csv")
      for i in range(len(df)):
          calc_conc_list.append( round(1e3 * pe / (TE.Na * TE.V_cyt),1) )
          calc_agg_list.append(df["Agg. coeff."][i])
          calc_speed_list.append(df["Speed"][i])
          calc_phitot_list.append(df["Phitot"][i])
          calc_phi_entire_list.append(pe)
          calc_Energy_list.append(df["Energy"][i])

    df_calc = pd.DataFrame()
    df_calc["conc"] = calc_conc_list
    df_calc["agg"] = calc_agg_list
    df_calc["speed"] = calc_speed_list
    df_calc["phitot"] = calc_phitot_list
    df_calc["phi_entire"] = calc_phi_entire_list
    df_calc["Energy"] = calc_Energy_list



    
    #plotting phitot deviation as heatmap
    calc_conc_list = []
    calc_agg_list = []
    calc_phitot_list = []
    calc_diffPhitot_list = []
    calc_Energy_list = []
    for pe in phi_entire_list:
      for agg in agg_list:
        print(pe, agg)
        p_calc = df_calc[df_calc["phi_entire"] == pe]
        k_calc = p_calc[p_calc["agg"] == agg]
        #p_sim = df_sim_heatmap[df_sim_heatmap["phi_entire"] == pe]
        #k_sim = p_sim[p_sim["agg"] == agg]
        calc_conc_list.append(list(k_calc["conc"])[0]) 
        calc_agg_list.append(list(k_calc["agg"])[0]) 
        calc_phitot_list.append(list(k_calc["phitot"])[0]) 
        calc_Energy_list.append(list(k_calc["Energy"])[0]) 
        
    df_calc_heatmap = pd.DataFrame()
    df_calc_heatmap["conc"] = calc_conc_list 
    df_calc_heatmap["agg"] = calc_agg_list 
    df_calc_heatmap["phitot"] = calc_phitot_list
    df_calc_heatmap["Energy"] = calc_Energy_list
    df_calc_heatmap.to_csv("phitot_calc.csv")
    res_calc_phitot = pivot_df(df_calc_heatmap, "agg", "conc", "phitot")
    ax41 = sns.heatmap(res_calc_phitot, ax = ax41, cmap = 'cividis')
    ax41.set_yticklabels(ax41.get_yticklabels(), rotation = 0)
    ax41.set_ylabel("$K_{agg}$")
    ax41.set_xlabel("Conc. $\mu M$")
    ax41.set_title("$\phi_{tot}^{calc}$")
    print("Calc phitot: ", res_calc_phitot)



"""
sim_conc_list = []
sim_agg_list = []
sim_phitot_list = []
sim_size_list = []
sim_phi_entire_list = []
sim_energy_list = []
for i in range(len(df_sim)):
    sim_conc = round(1e3 * df_sim["phi_entire"][i] / ( TE.Na * TE.V_cyt ), 1)
    sim_conc_list.append(sim_conc)
    sim_agg_list.append(df_sim["Aggregation"][i])
    sim_phitot_list.append(df_sim["phitot"][i])
    sim_size_list.append(df_sim["protSize"][i] * 1e6)
    sim_phi_entire_list.append(df_sim["phi_entire"][i])
    sim_energy_list.append(df_sim["energy"][i])
df_sim_heatmap = pd.DataFrame()
df_sim_heatmap["conc"] = sim_conc_list   
df_sim_heatmap["agg"] = sim_agg_list   
df_sim_heatmap["phitot"] = sim_phitot_list
df_sim_heatmap["size"] = sim_size_list   
df_sim_heatmap["phi_entire"] = sim_phi_entire_list
df_sim_heatmap["energy"] = sim_energy_list
print("Sim data")
print(df_sim_heatmap)


plot_45_sim = line_plots(["phi_entire", phi_entire_list], ["agg", select_kagg], df_sim_heatmap, "phitot")
plot_45_sim_energy = line_plots(["phi_entire", phi_entire_list], ["agg", select_kagg], df_sim_heatmap, "energy")
print("SIM ENERGY: ", plot_45_sim_energy)
wave_speed_phitot = list(pd.read_csv("../wave_sol/minData" + "_" + str(select_kagg) + ".csv")["phitot"])
wave_speed_energy = list(pd.read_csv("../wave_sol/minData" + "_" + str(select_kagg) + ".csv")["energy"])
print("Wave speed phitot: ", wave_speed_phitot)
conc_list = 1e3 * np.asarray(phi_entire_list) / (TE.Na * TE.V_cyt)
print("Plotting Phitot: ", plot_45_sim)    
#ax32.plot(conc_list, plot_45_sim, marker = "*", label = "Sim.")
#ax32.plot(conc_list, wave_speed_phitot, marker = "X", label = "$W_{s}$ analogy")
#ax32.plot(conc_list, landscape_min_phitot, "k", marker = "o", label = "$E_{min}$")
ax32.plot(conc_list, np.asarray(plot_45_sim_energy) / 1e3, marker = "+", label = "Sim.")
ax32.plot(conc_list, np.asarray(landscape_min_energy) / 1e3, "k", marker = "o", label = "$E_{min}$")
ax32.plot(conc_list, np.asarray(wave_speed_energy) / 1e3, marker = "X", label = "$W_{s}$ analysis")
ax32.set_xlabel("Conc $\mu M$")
ax32.set_ylabel("E ($10^{3} K_{B}T$)")
ax32.legend(frameon = False)
ax32.spines[['top', 'right']].set_visible(False)
ax32.set_xlim([0, 1])

#The following commented lines are plotting energies from min energy and wave speed min energy

"""    
minE, WaveE = GM.plot_WaveVsminEn() 
ax32.plot(minE[0], np.asarray(minE[1]) / 1e3, label = "$E_{min}$")
ax32.plot(WaveE[0], np.asarray(WaveE[1]) / 1e3, label = "$E_{W_{s}}$")
#ax32t = ax32.twinx()
error = (abs(np.asarray(minE[1]) - np.asarray(WaveE[1])) / abs(np.asarray(WaveE[1]))) * 100
#ax32t.plot(minE[0], error)
ax32.set_xlabel("Conc $\mu M$")
ax32.set_ylabel("$10^{3} K_{B}T$")
ax32.legend(frameon = False)
ax32.spines[['top', 'right']].set_visible(False)
"""

res = pivot_df(df_sim_heatmap, "agg", "conc", "phitot")
ax41 = sns.heatmap(res, ax = ax41, cmap = 'cividis', cbar_kws={'label': '$\phi_{tot}$'})
ax41.set_yticklabels(ax41.get_yticklabels(), rotation = 0)
ax41.set_ylabel("$K_{agg}$")
ax41.set_xlabel("Conc. $\mu M$")
    
print("SIZES: ", df_sim_heatmap["size"])
#Adding -1 for NaN values
for di in range(len(df_sim_heatmap["size"])):
      if np.isnan(df_sim_heatmap.iloc[di]["size"]) == True:
          df_sim_heatmap.iloc[di]["size"] = 0
print("SIZES: ", df_sim_heatmap["size"])
res = pivot_df(df_sim_heatmap, "agg", "conc", "size")
ax42 = sns.heatmap(res, ax = ax42, cmap = 'cividis', annot = False,  cbar_kws={'label': '2 X $r_{0}$'})
#ax42.set_ylabel("$K_{agg}$")
ax42.set_yticks([])
ax42.set_ylabel("")
ax42.set_xlabel("Conc. $\mu M$")
hatch_z = np.ma.masked_outside(res.values, -1, 0)
#hatch_z = np.ma.masked_invalid(res.values)
hatch_x = np.arange(len(res.columns) + 1)
hatch_y = np.arange(len(res.index) + 1)
ax42.pcolor(hatch_x, hatch_y, hatch_z, hatch = "//", alpha = 0)
plt.subplots_adjust(top = 0.95, bottom = 0.06, hspace = 0.4, wspace = 0.4)
plt.show()

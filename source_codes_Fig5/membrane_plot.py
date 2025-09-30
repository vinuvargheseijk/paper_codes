import numpy as np
import xml.etree.ElementTree as ET
import readXML
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
sns.set_context("poster")

Length = 50e-6
diffL = 20e-9

def extract_mem(trial_num, kcat, init_t, final_t, directory, axis):
    tag_string = "NstimL_10.0_kcat_enz_tiam_" + str(kcat) + "_spine_spacing_30.0_spacing_1.0_bfreq_20.0_coactive_1.0_mixed_True_trial_" + str(trial_num)
    cytosol = tag_string + "membrane.xml"
    chemTimes = pd.read_csv( directory + "/" + tag_string + "chemTimes.csv")["time"]
    tree1 = ET.parse(directory + "/" + cytosol)
    tree = [tree1]
    sum_an,max_an_y,an_y=readXML.plotXML(cytosol,tree)
    x_axis = np.arange(0, Length * 1e6, diffL * 1e6)
    x_heat = []
    value_heat = []
    time_heat = []
    time_init = chemTimes[chemTimes > init_t]
    time_final = chemTimes[chemTimes < final_t]
    time_init_index = list(time_init.index)[0]
    time_final_index = list(time_final.index)[-1]
    low_x = 0
    high_x = int(Length / diffL)
    print("Lower and Final time: ", time_init, time_init_index, time_final, time_final_index)
    for t in range(time_init_index, time_final_index):
        for x_in in range(low_x, high_x):
           time_heat.append(round(chemTimes.iloc[t], 1))
           x_heat.append(round(x_axis[x_in], 2))
           if an_y[t][x_in] < 1e-2:
              value_heat.append(np.nan)
           else:
              value_heat.append(an_y[t][x_in])

    df = pd.DataFrame()
    df["x"] = x_heat
    df["Time (s)"] = time_heat
    df["value"] = value_heat
    res = df.pivot(index = "Time (s)", columns = "x", values = "value")
    axis = sns.heatmap(res, ax = axis, cmap = 'cividis')
    axis.invert_yaxis()
    axis.set_xticks([])
    axis.set_xlabel("")
    ASx_start = int(10e-6/diffL)
    ASx_end = int(30e-6 / diffL)
    Sx_start = int(25e-6/diffL)
    Sx_end = int(10e-6 / diffL)
    rect1 = patches.Rectangle((ASx_start, 0), ASx_end, 20, linewidth=1, edgecolor='r', facecolor='orange')
    rect2 = patches.Rectangle((Sx_start, 0), Sx_end, 20, linewidth=1, edgecolor='r', facecolor='blue')
    axis.add_patch(rect1)
    axis.add_patch(rect2)


fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
extract_mem(0, 50.0, 0, 60, "./", ax1)
#extract_mem(5, 50.0, 0, 60, "./", ax3)
#extract_mem(10, 50.0, 0, 60, "./", ax5)
extract_mem(0, 200.0, 0, 60, "./", ax3)
#extract_mem(5, 200.0, 0, 60, "./", ax4)
#extract_mem(10, 200.0, 0, 60, "./", ax6)

plt.show()

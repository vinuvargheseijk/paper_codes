import matplotlib.pyplot as plt
import numpy as np
import dendShape
import pandas as pd
from mpl_toolkits import mplot3d

def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

fig = plt.figure()
ax = plt.axes(projection = '3d')

delay = 0.0
conc = 1.4e-3
spacing = 0.7
sstart = "CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "sstart.csv"
hmean = "CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "hmean.csv"
theta = "CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "theta.csv"
rp = "CytConc_" + str(conc) + "_ns_2_spine_create_delay_" + str(delay) + "_spine_spacing_" + str(spacing) + "rp.csv"
df_start = pd.read_csv("./" + sstart)
df_hmean = pd.read_csv("./" + hmean)
df_theta = pd.read_csv("./" + theta)
df_rp = pd.read_csv("./" + rp)

count = 0
for shape in range(0, len(df_start), 10):
   start_list = remove_nan(list(df_start.iloc[shape]))[1:]
   hmean_list = remove_nan(list(df_hmean.iloc[shape]))[1:]
   theta_list = remove_nan(list(df_theta.iloc[shape]))[1:]
   rp_list = remove_nan(list(df_rp.iloc[shape]))[1:]
   x, y = dendShape.get_shape(start_list, theta_list, hmean_list, rp_list)
   ax.plot3D(x, y, count * 0.2, color = 'g', alpha = max(y) / 0.25 )
   ax.set_xlim(10 - 0.6 * spacing, 10 + 0.6 * spacing)
   ax.set_ylim(0.0, max(y))
   count += 1

plt.show()



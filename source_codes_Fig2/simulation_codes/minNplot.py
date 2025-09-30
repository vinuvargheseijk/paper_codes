import pandas as pd
import matplotlib.pyplot as plt
import heatmap_minNplot as HM
import numpy as np



cyt_diff_list = HM.cyt_diff_list
phi_entire_list = HM.phi_entire_list



fig, (ax1, ax2) = plt.subplots(2)

for pe in phi_entire_list:
    df = pd.read_csv("minN" + str(pe) + ".csv")
    conc = round((pe / (HM.V_cyt * HM.Na)) * 1e3, 1)
    ax1.plot(np.asarray(df["cytD"]) * 1e12, df["minN"], marker = "*", label = "Conc: " + str(conc) + " $\mu M$" )
    ax1.legend()  
    #ax1.set_xlabel("Diffusion $\mu m^{2}/s$")
    ax1.set_ylabel("Minimum N")
    #ax1.set_title("Conc. $\mu M$  = " + str(conc))

    ax2.plot(np.asarray(df["cytD"]) * 1e12, df["minE"], marker = "*")
    ax2.set_xlabel("Diffusion $\mu m^{2}/s$")
    ax2.set_ylabel("Min. energy (KBT)")

plt.show()    

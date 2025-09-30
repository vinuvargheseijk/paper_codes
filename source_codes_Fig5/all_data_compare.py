import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
phisat = 1/50e-18
class data_extract:

  def data_save_lists(pe, directory, var_name):
    for ns in range(10):
      globals()[var_name + str(ns+1) + "_all"] = []
      globals()["phitot" + str(ns+1) + "_all"] = []
    var_field = pd.read_csv( directory + "/" + var_name + ".csv" )
    hmean = pd.read_csv( directory + "/" + "hmean" + ".csv" )
    theta = pd.read_csv( directory + "/" + "theta" + ".csv" )

    print("LENGTH of variable: ", len(var_field))
    for i in range(len(var_field)):  
      for j in range(1, len(var_field.iloc[i])):  
        globals()[var_name + str(j) + "_all"].append(var_field.iloc[i][j])
        area = 2 * np.pi * (1 / hmean.iloc[i][j])**2 * (1 - np.cos(theta.iloc[i][j]))
        phitot = var_field.iloc[i][j] * phisat * area
        globals()["phitot" + str(j) + "_all"].append(phitot)
    fin_phitot.append(globals()["phitot" + str(1) + "_all"][-1])    
    plt.plot(globals()["phitot" + str(1) + "_all"], label = "$\phi_{entire}$ = " + str(phi_entires[pe]))  
    plt.ylabel("$\phi_{tot}$")
    plt.legend()

phi_entires = [1000, 2000, 3000, 4000, 5000]
select_phitot = 10
filename = "phi.csv"
fin_phitot = []
phiE = []
for pe in range(len(phi_entires)):
    directory = "./phi_entire_" + str(phi_entires[pe]) + "/mns_" + str(select_phitot)
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == filename:
                print("File Found: ", file)
                data_extract.data_save_lists(pe, directory, "phi")
                phiE.append(phi_entires[pe])
                
plt.title("$\phi_{tot}$ = " + str(select_phitot))                
plt.show()
df = pd.DataFrame()
df["phi_entire"] = phiE
df["final_phitot"] = fin_phitot
df.to_csv("converged_values.csv")

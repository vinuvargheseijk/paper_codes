import pandas as pd
import numpy as np

def get_parameters(conc, spacing, spine_index):
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_25.0_spine_spacing_" + str(spacing) + "hmean.csv")
   print("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_25.0_spine_spacing_" + str(spacing) + "hmean.csv")
   print(df)
   for i in range(len(df)):
      if np.isnan(df.iloc[i]["1"]) == False:
         init_loc = i
         break
   df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_25.0_spine_spacing_" + str(spacing) + "shapeTimes.csv")
   print("Time: ", df.iloc[init_loc]["0"])
   params = ["hmean", "theta", "sstart", "rp"]
   for p in params:
      df = pd.read_csv("./CytConc_" + str(conc) + "_ns_2_spine_create_delay_25.0_spine_spacing_" + str(spacing) + p + ".csv")
      globals()[p] = df.iloc[init_loc][str(spine_index)]
      print(globals()[p])
   return globals()["hmean"], globals()["theta"], globals()["sstart"], globals()["rp"]   

hmean, theta, sstart, rp = get_parameters(0.001, 1.5, 1)
 

import pandas as pd
import numpy as np
import dendShape
import writeXML
import matplotlib.pyplot as plt
import totalEnergy
#import chemModel as CM
import readXML as RX
import xml.etree.ElementTree as ET
import totalEnergy as TE
import tag
import seaborn as sns
from scipy import optimize

def normalize_exp_sim(data):
    return np.asarray(data) / max(data)

def shift_time_exp(time_exp, shift):
    time_exp = np.asarray(time_exp) * 1e-3 + shift
    return time_exp

sns.set_context("poster")
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
       

V_voxel = np.pi * 0.5 * (TE.dendDia)**2 * TE.diffL        
dt = 0.001
plot_interval = 1
t_axis = np.linspace(0, tag.runtime, int(tag.runtime/dt))
df_exp = pd.read_csv("./Calcium_Helmchen.csv")
tag_string = tag.tag_string
CaXML = tag_string + "Ca.xml"
RacXML = tag_string + "Rac.xml"
cytosolXML = tag_string + "cytosol.xml"
thrXML = tag_string + "thr286.xml"
check_files = [CaXML, RacXML, cytosolXML, thrXML]
names = ["Ca", "Rac", "cytosol", "thr286"]
mol_name = 0
norm_exp = normalize_exp_sim(df_exp["conc"])
time_exp = df_exp["time"]
shift_time = 0.475 #simulation pulse comes at 0.5 s
time_exp = shift_time_exp(time_exp, shift_time)
ax2.plot(np.asarray(time_exp), norm_exp, label = "Exp.")
for c in check_files:
     filename = c
     tree1 = ET.parse(filename)
     tree = [tree1]
     sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
     if names[mol_name] == "Rac":
         #ax1.plot( t_axis, normalize_exp_sim(max_an_y), 'r', label = names[mol_name] )
         ax1.plot( t_axis, 1e3 * np.asarray(max_an_y) / (V_voxel * TE.Na) , 'r', label = names[mol_name] )
     if names[mol_name] == "Ca":
         #ax2.plot( t_axis, normalize_exp_sim(np.asarray(max_an_y) - min(max_an_y)), 'b', label = names[mol_name] )
         #ax2.plot( t_axis, 1e3 *  np.asarray(max_an_y) / (V_voxel * TE.Na), 'b', label = names[mol_name] ) #Actual conc.
         norm_Ca = normalize_exp_sim(np.asarray(max_an_y)- min(max_an_y))
         ax2.plot( t_axis, norm_Ca, 'b', label = names[mol_name] ) #Actual conc.
     if names[mol_name] == "cytosol":
         #ax1.plot( t_axis, normalize_exp_sim(max_an_y), 'b', label = names[mol_name] )
         ax3.plot( t_axis, np.asarray(sum_an), 'g', label = names[mol_name] )
     if names[mol_name] == "thr286":
         #ax1.plot( t_axis, normalize_exp_sim(max_an_y), 'b', label = names[mol_name] )
         ax3.plot( t_axis, np.asarray(sum_an) , 'c', label = names[mol_name] )
     mol_name = mol_name + 1
     ax2.set_title( "Ca sim. conc") 
     ax2.set_xlabel("Time (s)")
     ax2.legend()
     ax3.set_title( "Cytosol, Membrane #") 
     ax3.legend()
     ax1.legend()


plt.show()     

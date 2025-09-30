import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import dendShape
import seaborn as sns
sns.set_context("poster")
#sns.set(font_scale = 0.5)

Length = 20e-6
Rc = 0.5e-6
V_cyt = np.pi * 0.5e-6**2 * Length
Na = 6.023e23
phisat =  1/50e-18

figH = plt.figure(figsize = (14, 48))
ax11H = plt.subplot2grid((4, 2), (0, 0), colspan = 1)
ax12H = plt.subplot2grid((4, 2), (0, 1), colspan = 1)
ax21H = plt.subplot2grid((4, 2), (1, 0), colspan = 1)
ax22H = plt.subplot2grid((4, 2), (1, 1), rowspan = 1)
ax31H = plt.subplot2grid((4, 2), (2, 0), rowspan = 1)
ax32H = plt.subplot2grid((4, 2), (2, 1), rowspan = 1)
ax41H = plt.subplot2grid((4, 2), (3, 0), rowspan = 1)
ax42H = plt.subplot2grid((4, 2), (3, 1), rowspan = 1)
#ax51H = plt.subplot2grid((5, 2), (4, 0), rowspan = 1)
#ax52H = plt.subplot2grid((5, 2), (4, 1), rowspan = 1)
#ax61H = plt.subplot2grid((6, 2), (5, 0), rowspan = 1)
#ax62H = plt.subplot2grid((6, 2), (5, 1), rowspan = 1)

ax11H.text(-0.1, 1.1, 'A', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax11H.transAxes)


ax11H.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
ax12H.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
ax11H.set_xticks([])
ax12H.set_xticks([])
ax11H.set_yticks([])
ax12H.set_yticks([])


ax12H.text(-0.1, 1.1, 'B', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax12H.transAxes)

ax21H.text(-0.1, 1.1, 'C', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax21H.transAxes)

ax22H.text(-0.1, 1.1, 'D', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax22H.transAxes)

ax31H.text(-0.1, 1.1, 'E', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax31H.transAxes)

ax32H.text(-0.1, 1.1, 'F', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax32H.transAxes)

ax41H.text(-0.1, 1.1, 'G', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax41H.transAxes)

ax42H.text(-0.1, 1.1, 'H', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax42H.transAxes)
#ax51H.text(-0.1, 1.1, 'I',
#        horizontalalignment='left',
#        verticalalignment='bottom',
#        transform=ax51H.transAxes)
#ax52H.text(-0.1, 1.1, 'J',
#        horizontalalignment='left',
#        verticalalignment='bottom',
#        transform=ax52H.transAxes)
"""
ax61H.text(-0.1, 1.1, 'K',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax61H.transAxes)
ax62H.text(-0.1, 1.1, 'L',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax62H.transAxes)
"""

shape_colors = ["b", "g", "k", "m", "c"]
bar_colors = ['tab:red', 'tab:blue', 'tab:red', 'tab:orange']
def low_phi_absolute():
   #phi_entire = 1800.0
   phi_entire = 1200.0
   k_agg1 = -60.0
   return phi_entire, k_agg1

def high_phi_absolute():
   phi_entire = 6000.0
   k_agg1 = -60.0
   return phi_entire, k_agg1

phi_entire, k_agg1 = low_phi_absolute()

df_branches = pd.read_csv("./mu-4.5_probability_" + str(k_agg1) + "_" + str(phi_entire) + "_5.0nm.csv")
low_conc = round(1e3 * phi_entire / (Na * V_cyt), 1)
ax31H.plot(df_branches["phitot"], np.asarray(df_branches["lowEnergy"]) / 1e3, 'k', linestyle = "--")
ax31H.plot(df_branches["phitot"], np.asarray(df_branches["highEnergy"]) / 1e3, 'k', linestyle = "-")
LowEnergy_coord = list(df_branches["lowEnergy"])
minLowEnergy_coord_index = LowEnergy_coord.index(np.nanmin(LowEnergy_coord))
ax31H.text(df_branches.iloc[minLowEnergy_coord_index]["phitot"], np.nanmin(LowEnergy_coord)/1e3 - 0.2, "$E_{Shallow}$")
ax31H.plot(df_branches.iloc[minLowEnergy_coord_index]["phitot"], np.nanmin(LowEnergy_coord)/1e3, marker = "X", markerfacecolor = "k", markeredgecolor = "k")
HighEnergy_coord = list(df_branches["highEnergy"])
minHighEnergy_coord_index = HighEnergy_coord.index(np.nanmin(HighEnergy_coord))
ax31H.text(df_branches.iloc[minHighEnergy_coord_index]["phitot"]+200, np.nanmin(HighEnergy_coord)/1e3 - 0.1, "$E_{Sharp}$")
ax31H.plot(df_branches.iloc[minHighEnergy_coord_index]["phitot"], np.nanmin(HighEnergy_coord)/1e3, marker = "o", markerfacecolor = "k", markeredgecolor = "k")
low_phi_phitot = df_branches.iloc[minLowEnergy_coord_index]["phitot"]
high_phi_phitot = df_branches.iloc[minHighEnergy_coord_index]["phitot"]
print("HCoords: ", df_branches.iloc[minHighEnergy_coord_index]["phitot"], np.nanmin(HighEnergy_coord)/1e3)
print("LCoords: ", df_branches.iloc[minLowEnergy_coord_index]["phitot"], np.nanmin(LowEnergy_coord)/1e3)
ax31H.set_ylim(-2, 0)
ax31H.set_ylabel("$10^{3} * K_{B}T$")
ax31H.set_title("Total conc. = " + str(low_conc) + "$\mu M$")
ax31H.set_xlabel("$\phi_{tot}$")
ax31H.spines[['right', 'top']].set_visible(False)
ax31H.legend(frameon = False)
ax31H.set_xlim([0, max(df_branches["phitot"])])


phi_entire, k_agg1 = high_phi_absolute()
###Reading phitots in the geom_param_file
df = pd.read_csv("./" + str(phi_entire) + "_" + str(k_agg1) + "geom_param.csv")

##Plotting branches for high phi
df_branches = pd.read_csv("./mu-4.5_probability_" + str(k_agg1) + "_" + str(phi_entire) + "_5.0nm.csv")

high_conc = round(1e3 * phi_entire / (Na * V_cyt), 1)
ax32H.yaxis.set_label_position("right")
ax32H.plot(df_branches["phitot"], np.asarray(df_branches["lowEnergy"]) / 1e3, 'k', linestyle = "--")
ax32H.plot(df_branches["phitot"], np.asarray(df_branches["highEnergy"]) / 1e3, 'k', linestyle = "-")
ax32H.set_ylim(-20, 0)
#ax32H.set_ylabel("$K_{B}T$")
LowEnergy_coord = list(df_branches["lowEnergy"])
minLowEnergy_coord_index = LowEnergy_coord.index(np.nanmin(LowEnergy_coord))
ax32H.text(df_branches.iloc[minLowEnergy_coord_index]["phitot"], np.nanmin(LowEnergy_coord)/1e3 - 0.2, "$E_{Shallow}$")
ax32H.plot(df_branches.iloc[minLowEnergy_coord_index]["phitot"], np.nanmin(LowEnergy_coord)/1e3, marker = "X", markerfacecolor = "k", markeredgecolor = "k")
ax32H.set_xlabel("$\phi_{tot}$")
ax32H.spines[['right', 'top']].set_visible(False)
ax32H.set_title("Total conc. = " + str(high_conc) + "$\mu M$")
ax32H.legend(frameon = False)
HighEnergy_coord = list(df_branches["highEnergy"])
minHighEnergy_coord_index = HighEnergy_coord.index(np.nanmin(HighEnergy_coord))
ax32H.text(df_branches.iloc[minHighEnergy_coord_index]["phitot"] - 200, np.nanmin(HighEnergy_coord)/1e3 - 1.25, "$E_{Sharp}$")
ax32H.plot(df_branches.iloc[minHighEnergy_coord_index]["phitot"], np.nanmin(HighEnergy_coord)/1e3, marker = "o", markerfacecolor = "k", markeredgecolor = "k")
ax32H.set_xlim([0, max(df_branches["phitot"])])


def line_plot_energy():
   #phi_entires = [1000.0, 1400.0, 1800.0, 2000.0, 3000.0, 4000.0, 6000.0]
   phi_entires = [500.0, 1200.0, 1400.0, 1800.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 8000.0]
   k_agg1 = -60.0
   LCON = []
   LEN = []
   LAG = []
   LME = []
   LMPot = [] 
   LCPot = []
   LRm = []
   LCtotalPot = []
   LTotal = []
   HCON = []
   HEN = []
   HAG = []
   HME = [] 
   HMPot = []
   HCPot = []
   HRm = []
   HCtotalPot = []
   HTotal = []
   Hsort = []
   Lsort = []
   for pe in phi_entires:
       print("PHI ENTIRE: ", pe)
       df = pd.read_csv("./" + str(pe) + "_" + str(k_agg1) + "geom_param.csv")
       minLowEnergy = df[df["lowTotal"] == np.nanmin(list(df["lowTotal"])[0:-2])]
       LME.append(minLowEnergy["LME"]) 
       LEN.append(minLowEnergy["LEN"]) 
       LAG.append(minLowEnergy["LAG"]) 
       LCON.append(minLowEnergy["LCON"])
       LMPot.append(minLowEnergy["LMPot"])
       LCE_per_molecule = minLowEnergy["LCE"] / minLowEnergy["phitot"]
       LCPot.append(LCE_per_molecule)
       LRm.append(minLowEnergy["LRm"])
       LCtotalPot.append(LCE_per_molecule + minLowEnergy["LMPot"])
       LTotal.append(minLowEnergy["lowTotal"])    
       
       lowArea_t = 2 * np.pi * (minLowEnergy["LRm"])**2 * (1 - np.cos(minLowEnergy["LTheta"]))
       lowR0 = 2 * (minLowEnergy["LRm"] + minLowEnergy["LRp"]) * np.sin(minLowEnergy["LTheta"])
       lowArea_c = 2 * np.pi * Rc * Length - np.pi * lowR0**2
       low_phi_t = minLowEnergy["phitot"] / (phisat * lowArea_t)
       low_phi_c = (pe - float(minLowEnergy["phitot"])) / (phisat * lowArea_c)
       Lsort.append(low_phi_t/ low_phi_c)
       
       minHighEnergy = df[df["highTotal"] == np.nanmin(list(df["highTotal"])[0:-2])]
       HME.append(minHighEnergy["HME"]) 
       HEN.append(minHighEnergy["HEN"]) 
       HAG.append(minHighEnergy["HAG"]) 
       HCON.append(minHighEnergy["HCON"])
       HMPot.append(minHighEnergy["HMPot"])
       HRm.append(minLowEnergy["HRm"])
       HTotal.append(minLowEnergy["highTotal"])    
       HCE_per_molecule = np.asarray(minHighEnergy["HCE"]) / np.asarray(minHighEnergy["phitot"])
       HCtotalPot.append(HCE_per_molecule + minHighEnergy["HMPot"])
       HCPot.append(HCE_per_molecule)
       highArea_t = 2 * np.pi * (minHighEnergy["HRm"])**2 * (1 - np.cos(minHighEnergy["HTheta"]))
       highR0 = 2 * (minHighEnergy["HRm"] + minHighEnergy["HRp"]) * np.sin(minHighEnergy["HTheta"])
       highArea_c = 2 * np.pi * Rc * Length - np.pi * highR0**2
       high_phi_t = minHighEnergy["phitot"] / (phisat * highArea_t)
       high_phi_c = (pe - float(minHighEnergy["phitot"])) / (phisat * highArea_c)
       print("Surface fractions: ", high_phi_t, high_phi_c)
       Hsort.append(high_phi_t/ high_phi_c)
   print("Mech potential: ", LMPot)
   print(HMPot)
   conc_list = 1e3 * np.asarray(phi_entires) / (Na * V_cyt)    
   #ax51H.plot(conc_list, Hsort, "k", marker = "o", markerfacecolor = "k", markeredgecolor = "k", label = "High $\phi$")
   #ax51H.plot(conc_list, Lsort, "k", marker = "o", markerfacecolor = "None", markeredgecolor = "k", label = "Low $\phi$")
   #ax51H.set_ylabel("Sorting")
   #ax51H.legend(frameon = False)
   #ax52H.plot(conc_list, np.asarray(HCPot) + np.asarray(HMPot), "k", label = "$\mu_{total}$")
   #ax52H.set_ylabel("$K_{B}T$")
   #ax52H.set_xlabel("Conc $\mu M$")
   #ax51H.set_xlabel("Conc $\mu M$")
   #ax52H.legend()
   ax21H.plot(conc_list, np.asarray(LME) / 1e3, "k", label = "Mismatch", marker = "o")    
   ax21H.plot(conc_list, np.asarray(LAG) / 1e3, "k", label = "Aggregation", marker = "x")    
   ax21H.plot(conc_list, np.asarray(LEN) / 1e3, "k", label = "Entropy", marker = "*")    
   ax21H.plot(conc_list, np.asarray(LCON) / 1e3, "k", label = "Constraint", marker = "v")    
   ax21H.plot(conc_list, [0] * len(conc_list), 'k', linestyle = "--")
   print(LME)
   print(HME)
   ax21H.set_xlabel("Conc $\mu M$")
   ax21H.set_ylabel("$10^{3} * K_{B}T$")
   
   print(HME)
   ax22H.plot(conc_list, np.asarray(HME) / 1e3, "k", label = "Mismatch", marker = "o")    
   ax22H.plot(conc_list, np.asarray(HAG) / 1e3, "k", label = "Aggregation", marker = "x")    
   ax22H.plot(conc_list, np.asarray(HEN) / 1e3, "k", label = "Entropy", marker = "*")    
   ax22H.plot(conc_list, np.asarray(HCON) / 1e3, "k", label = "Constraint", marker = "v")    
   #ax42H.plot(conc_list, np.asarray(HCPot), "r", label = "Mpotential")    
   #ax42H.plot(conc_list, np.asarray(HCtotalPot), "k", label = "Mpotential")    
   ax22H.plot(conc_list, [0] * len(conc_list), 'k', linestyle = "--")
   ax22H.set_xlabel("Conc $\mu M$")
   #ax22H.set_ylabel("$10^{3} * K_{B}T$")
   ax21H.spines[['right', 'top']].set_visible(False)
   ax22H.spines[['right', 'top']].set_visible(False)
   
   #ax21H.legend()
line_plot_energy()

phi_entire, k_agg1 = low_phi_absolute()
df = pd.read_csv("./" + str(phi_entire) + "_" + str(k_agg1) + "geom_param.csv")

#The following commented code is for plotting the shape variation with respect to phitot change
"""
shape_phitots = [ df.iloc[i]["phitot"] for i in range(0, len(list(df["phitot"])), 5) ]
shape_indices = [ i for i in range(0, len(list(df["phitot"])), 5) ]

for sptot in range(len(shape_phitots) - 1):
    ax21H.arrow(shape_phitots[sptot], -1, 0, 0.1, width = 12.5, head_length = 0.1, facecolor = "k", alpha = shape_phitots[sptot] / max(shape_phitots))
"""

def translate_center_shape(x):
  x = list(x)
  mid_value = 0.5 * (x[-1] - x[0])
  print(x[-1])
  print(x[0])
  print(mid_value)
  mid_point = np.abs(np.asarray(x) - mid_value).argmin()
  print(mid_point)
  x = np.asarray(x) - x[mid_point]
  return x, x[mid_point]


#Finding index in the geom files. The low of low phi and high are stored as -2 and -1 indices. Lets check this
print("Low phi phitot: ", low_phi_phitot)
print("Low phi phitot from geom file: ", df.iloc[-2]["phitot"])
print(df.iloc[-2]["LRm"])
x, y = dendShape.get_shape([0e-6], [df.iloc[-2]["LTheta"]], [1 / df.iloc[-2]["LRm"]], [df.iloc[-2]["LRp"]])
x, mid_point = translate_center_shape(x)
ax41H.plot(x,np.asarray(y) * 1e3, color = "k", linestyle = "--", label = "$E_{Shallow}$")

print("High phi phitot: ", high_phi_phitot)
print("High phi phitot from geom file: ", df.iloc[-1]["phitot"])
print(df.iloc[-1]["HRm"])
x, y = dendShape.get_shape([0e-6], [df.iloc[-1]["HTheta"]], [1 / df.iloc[-1]["HRm"]], [df.iloc[-1]["HRp"]])
x, mid_point = translate_center_shape(x)
ax41H.plot(x,np.asarray(y) * 1e3, color = "k", label = "$E_{Sharp}$")
ax41H.set_ylabel("Height nm")
ax41H.legend(frameon = False)
ax41H.spines[['right', 'top']].set_visible(False)
#The following commented code is for plotting the shape variation with respect to phitot change
"""
for i in range(len(shape_indices)):
  index = shape_indices[i]
  #rm = df.iloc[index]["LRm"] 
  #theta = df.iloc[index]["LTheta"] 
  #rp = df.iloc[index]["LRp"]
  #r0 = (rm + rp) * np.sin(theta)
  x, y = dendShape.get_shape([0e-6], [df.iloc[index]["LTheta"]], [1 / df.iloc[index]["LRm"]], [df.iloc[index]["LRp"]])
  x, mid_point = translate_center_shape(x)
  #ax41H.plot(x,np.asarray(y), color = "k", alpha = shape_phitots[i] / max(shape_phitots) )
  x, y = dendShape.get_shape([0e-6], [df.iloc[index]["HTheta"]], [1 / df.iloc[index]["HRm"]], [df.iloc[index]["HRp"]])
  x, mid_point = translate_center_shape(x)
  #ax51H.plot(x,np.asarray(y), color = shape_colors[i])
  ax51H.plot(x,np.asarray(y), color = "k", alpha = shape_phitots[i] / max(shape_phitots) )

#ax4.set_yticks([0])
ax51H.set_xlabel("x $\mu m$")
ax51H.set_ylabel("High $\phi$\nh $\mu m$")
"""
#ax51H.spines[['right', 'top']].set_visible(False)
#ax51H.spines[['right', 'top']].set_visible(False)

phi_entire, k_agg1 = high_phi_absolute()
df = pd.read_csv("./" + str(phi_entire) + "_" + str(k_agg1) + "geom_param.csv")

print(df.iloc[-1]["HRm"])
x, y = dendShape.get_shape([0e-6], [df.iloc[-1]["HTheta"]], [1 / df.iloc[-1]["HRm"]], [df.iloc[-1]["HRp"]])
x, mid_point = translate_center_shape(x)
ax42H.plot(x,np.asarray(y) * 1e3, color = "k", label = "$E_{Sharp}$")
ax42H.set_ylabel("Height nm")
ax42H.legend(frameon = False)
ax42H.spines[['right', 'top']].set_visible(False)


shape_phitots = [ df.iloc[i]["phitot"] for i in range(2, len(list(df["phitot"])), 5) ]
shape_indices = [ i for i in range(2, len(list(df["phitot"])), 5) ]

#The following commented code is for plotting the shape variation with respect to phitot change
"""
for sptot in range(len(shape_phitots)):
    ax22H.arrow(shape_phitots[sptot], -10, 0, 1, width = 50, head_length = 1, facecolor = "k", alpha = shape_phitots[sptot] / max(shape_phitots))


for i in range(len(shape_indices)):
  index = shape_indices[i]
  #rm = df.iloc[index]["LRm"] 
  #theta = df.iloc[index]["LTheta"] 
  #rp = df.iloc[index]["LRp"]
  #r0 = (rm + rp) * np.sin(theta)
  x, y = dendShape.get_shape([0e-6], [df.iloc[index]["LTheta"]], [1 / df.iloc[index]["LRm"]], [df.iloc[index]["LRp"]])
  x, mid_point = translate_center_shape(x)
  #ax42H.plot(x,np.asarray(y) * 1e3, color = shape_colors[i])
  #ax42H.plot(x,np.asarray(y), color = "k", alpha = shape_phitots[i] / max(shape_phitots) )
  x, y = dendShape.get_shape([0e-6], [df.iloc[index]["HTheta"]], [1 / df.iloc[index]["HRm"]], [df.iloc[index]["HRp"]])
  x, mid_point = translate_center_shape(x)
  #ax52H.plot(x,np.asarray(y))
  #ax52H.plot(x,np.asarray(y), color = shape_colors[i])
  ax52H.plot(x,np.asarray(y), color = "k", alpha = shape_phitots[i] / max(shape_phitots) )
"""
ax42H.spines[['right', 'top']].set_visible(False)
ax42H.spines[['right', 'top']].set_visible(False)
#ax52H.spines[['right', 'top']].set_visible(False)
#ax52H.spines[['right', 'top']].set_visible(False)

ax41H.set_xlabel("x $\mu m$")
ax42H.set_xlabel("x $\mu m$")

#circle00 = plt.Circle((0.0,-0.5), 0.5, color = 'grey')
#ax41H.add_patch(circle00)
#circle01 = plt.Circle((0.0,-0.5), 0.5, color = 'grey')
#ax42H.add_patch(circle01)

#circle1 = plt.Circle((0.0,-0.5), 0.5, color = 'grey')
#ax51H.add_patch(circle1)
#circle2 = plt.Circle((0.0,-0.5), 0.5, color = 'grey')
#ax52H.add_patch(circle2)
#ax52H.set_xlim([-0.5, 0.5])
#ax51H.set_xlim([-0.5, 0.5])
#ax52H.set_ylim([-0.05, 0.15])
#ax51H.set_ylim([-0.05, 0.15])
figH.subplots_adjust(
    hspace=0.4,
    wspace=0.4
)

print(df)
plt.show()

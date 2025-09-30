import pandas as pd
import numpy as np
import writeXML
import matplotlib.pyplot as plt
import readXML as RX
import xml.etree.ElementTree as ET
import seaborn as sns
from scipy import optimize
import energy_process as EP
import argparse
import matplotlib.patches as patches
sns.set_context("poster")

parser = argparse.ArgumentParser()
parser.add_argument("--mu0", type = float)
args = parser.parse_args()

KBT = EP.KBT
Na = EP.Na
RT = EP.RT
RT_kcal_mol = EP.RT_kcal_mol #kcal/mol
k = EP.k
khat = EP.khat
kentropy = KBT
Cp = EP.Cp
phisat = EP.phisat
phi_entire_list = []
k_agg_list = []
global state_list, energy_list, dE_list, phi_list, diff_phitot_list, phitot_list
state_list = []
sum_list = []
energy_list = []
dE_list = []
phi_list = []
diff_phitot_list = []
phitot_list = []
Ca_list = []
Rac_list = []
cytosol_list = []
mu0_exp = -4.0
scalF = 1.0
Kb = 40.0
dt = 0.24
cyt_diffConst = 0.5e-12
Length = 10e-6
Rc = 0.5e-6
Na = 6.023e23
V_cyt = np.pi * Rc**2 * Length 
V_voxel = np.pi * 0.5 * (1e-6)**2 * 20e-9        
CaTau = 0.08

bg_int_list = [1.0]
nb_list = [12]
num_bursts = 1
#@nb_list = [50, 75, 100]
bg_freq_list = [1.0]
freq_ratio_list = [1.0]
freq_list = [0.5, 1.0, 2.0, 4.0, 8.0]
#freq_list = [4.0, 8.0]
const_param1 = bg_freq_list[0]
const_param2 = freq_ratio_list[0]
const_param3 = bg_int_list[0]
source_list = [8000.0, 16000.0]
k_agg = -45.0
mu0_exp = -4.5
source = source_list[0]
phi_entire = 6000.0
y_pos_list= [0.5e-6, 1.0e-6, 1.5e-6, 2.0e-6]
y_pos = y_pos_list[0]
ypos_ratio = 0.0
spacing = 0.75e-6
var_param1 = y_pos_list
var_param2 = source_list
num_spines = 1

hmp_r1 = ["y"] #Check if const and var params are correct
hmp_r2 = ["#glut."]
ECD = 0.3e-12
Source = 8000

def glutamate_profile(height, time):
   elec_dx = 0.1e-6
   deltaY = elec_dx
   deltaZ = height * 1e-6
   conc_factor = 1 / (Na * height * 1e-6)
   factor = (elec_dx * deltaY * Na * deltaZ)  * conc_factor
   return factor * ( Source/ (4 * np.pi * time * ECD) ) * np.exp( -1 * ( (height * 1e-6)**2 ) / (4 * ECD * time) )

def line_plots(v_parameters, c_parameters, df, plot_name):
    print(df[c_parameters[0]][0])
    plot_list = []
    for i in range(len(df)):
        print(df[c_parameters[0]][i], c_parameters[1])
        if df[c_parameters[0]][i] == c_parameters[1]:
            print("Found: ", df[c_parameters[0]][i])
            for j in v_parameters[1]:
               print(j) 
               if df[v_parameters[0]][i] == j:
                   print("Found", df[v_parameters[0]][i])
                   print(df[plot_name][i])
                   plot_list.append(df[plot_name][i])
    return plot_list            

def find_max_phitot(phi, hmean, theta):
    area = 2 * np.pi * (1/np.asarray(hmean))**2 * (1 - np.cos(np.asarray(theta)))
    phitot = list( np.asarray(phi) * phisat * area )
    maxid = phitot.index(max(phitot))
    print(maxid)
    return maxid

def spine_df(tag_string): 
       try:
          print("reading file") 
          phi = pd.read_csv(tag_string + "phi.csv")
          theta = pd.read_csv(tag_string + "theta.csv")
          saddle_start = pd.read_csv(tag_string + "sstart.csv")
          hmean = pd.read_csv(tag_string + "hmean.csv")
          rp = pd.read_csv(tag_string + "rp.csv")
          mConc = pd.read_csv(tag_string + "mConc.csv")
          cConc = pd.read_csv(tag_string + "cConc.csv")
          total_energy = 0
          total_phitot = 0
          state = ""
          sum_t = 0
          phitot_ns = []
          #counting spines
          sp_count = 0
          for ns in range(1, 10):
              try:
                  phi.iloc[-1][ns]
                  sp_count = sp_count + 1
              except:
                  print("Num spines: ", sp_count)
          for ns in range(1, sp_count+1):
             try:  
               in_index = find_max_phitot(list(phi.iloc[:,ns]), list(hmean.iloc[:,ns]), list(theta.iloc[:,ns]))
               conc = cConc.iloc[in_index][ns]
               rm_last = 1 / hmean.iloc[in_index][ns]
               theta_last = theta.iloc[in_index][ns]
               init_area = EP.calc_dome_area(1/hmean.iloc[in_index][ns], theta.iloc[in_index][ns])
               area = EP.calc_dome_area(rm_last, theta_last)
               phi_last = phi.iloc[in_index][ns]
               init_phitot = phi.iloc[in_index][ns] * phisat * init_area
               phitot = phi_last * phisat * area
               phitot_ns.append(phitot)
               init_energy = EP.total_energy([(1/hmean.iloc[in_index][ns]) * 1e6 * theta.iloc[in_index][ns], theta.iloc[in_index][ns], rp.iloc[in_index][ns] * 1e6], init_phitot, cConc.iloc[in_index][ns], k_agg, mu0_exp)  #Check for input arguments
               energy = EP.total_energy([rm_last * 1e6 * theta_last, theta_last, rp.iloc[in_index][ns] * 1e6], phitot, conc, k_agg, mu0_exp) #Check for input arguments
               if np.isnan(energy) == False:
                 total_energy = total_energy + energy
                 total_phitot = total_phitot + phitot
                 
               if phi.iloc[in_index][ns] <= 0.1:
                   state = state + "0"
                   sum_t = sum_t + 0
               else:
                   state = state + "1"
                   sum_t = sum_t + 1
             except:
                   print("Number of spines: " + str(ns))
          energy_list.append(total_energy)
          state_list.append(state)
          phitot_list.append(total_phitot)
          try:
             diff_phitot_list.append(phitot_ns[0] - phitot_ns[1])
          except:
             diff_phitot_list.append(np.nan)  
          sum_list.append(sum_t)
          phi_list.append(phi.iloc[in_index][1])
       except Exception as excpt:
           type_excpt = type(excpt).__name__
           print(type_excpt)
           print(tag_string)
           if type_excpt != "IndexError":
             energy_list.append(np.nan)
             state_list.append(np.nan)
             diff_phitot_list.append(np.nan)
             phitot_list.append(np.nan)
             phi_list.append(np.nan)
             sum_list.append(np.nan)
       print(len(energy_list))
       print(len(diff_phitot_list))
       print(len(phi_entire_list))


for vp1 in var_param1:
    for vp2 in var_param2:
       tag_string = "Source_" + str(vp2) + "_nb_" + str(num_bursts) + "_ns_" + str(num_spines) + "_ypos_" + str(vp1) + "_spacing_" + str(spacing) + "_freq_" + str(1.0) + "_CaTau_" + str(CaTau) + "_phi_entire_" + str(phi_entire)
       CaXML = tag_string + "Ca.xml"
       RacXML = tag_string + "Rac.xml"
       cytosolXML = tag_string + "cytosol.xml"
       check_files = [CaXML, RacXML, cytosolXML]
       print("TAG STRING: ", tag_string)
       hmp_r1.append(vp1)
       hmp_r2.append(vp2)
       for c in check_files:
            filename = c
            try:
                tree1 = ET.parse(filename)
            except:
                if c == CaXML:
                   Ca_list.append(np.nan)
                if c == RacXML:
                   Rac_list.append(np.nan)
                if c == cytosolXML:
                   cytosol_list.append(np.nan)
                break
            tree = [tree1]
            sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
            try:
              if c == CaXML:
                idy = sum_an.index(max(sum_an))
                print(idy)
                Ca_list.append(1e3 * np.sum(an_y[idy]) / (Na * V_cyt))
              if c == RacXML:
                idy = sum_an.index(max(sum_an))
                Rac_list.append(1e3 * np.sum(an_y[idy]) / (Na * V_cyt))
              if c == cytosolXML:
                idy = sum_an.index(max(sum_an))
                cytosol_list.append(np.sum(an_y[idy]))           
            except:     
              if c == CaXML:
                Ca_list.append(np.nan)
              if c == RacXML:
                Rac_list.append(np.nan)
              if c == cytosolXML:
                cytosol_list.append(np.nan)           
            
df = pd.DataFrame()
df[hmp_r1[0]] = hmp_r1[1:]
df[hmp_r2[0]] = hmp_r2[1:]
print(df)
df["Ca"] = Ca_list
print(df)

"""
try:
   pd.read_csv("allParams.csv")
   df.to_csv("allParams.csv", mode = 'a', header = False, index = False)
except:
   df.to_csv("allParams.csv", index = False)
"""

#Plotting the calcium Sim Vs Exp.
figH = plt.figure(figsize = (10, 18))
ax11H = plt.subplot2grid((3, 2), (0, 0), colspan = 1)
ax12H = plt.subplot2grid((3, 2), (0, 1), colspan = 1)
ax21H = plt.subplot2grid((3, 2), (1, 0), colspan = 1)
ax22H = plt.subplot2grid((3, 2), (1, 1), rowspan = 1)
ax31H = plt.subplot2grid((3, 2), (2, 0), colspan = 1)
ax32H = plt.subplot2grid((3, 2), (2, 1), rowspan = 1)

ax11H.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax12H.spines[['top', 'right']].set_visible(False)
ax21H.spines[['top', 'right']].set_visible(False)
ax22H.spines[['top', 'right']].set_visible(False)
ax31H.spines[['top', 'right']].set_visible(False)
ax32H.spines[['top', 'right']].set_visible(False)
ax11H.set_xticks([])
ax11H.set_yticks([])

ax11H.text(-0.1, 1.1, 'A', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax11H.transAxes)

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


def RMSD_calc(exp_time, sim_time, data_exp, data_sim):
   tol = 1e-1
   indices_exp = []
   indices_sim = []
   time_exp = []
   time_sim = []
   sum_sq = 0
   for i in range(len(exp_time)):
       for j in range(len(sim_time)):
           if abs(sim_time[j] - exp_time[i]) < tol:
              min_diff = sim_time[j]
              index_min = j
              indices_sim.append(j)
              indices_exp.append(i)
              time_exp.append(exp_time[i])
              time_sim.append(sim_time[j])
              break
       sum_sq = sum_sq + (abs(data_sim[j] - data_exp[i]))**2
   print("TIMES (exp) for RMSD: ", time_exp)
   print("TIMES (sim) for RMSD: ", time_sim)
   return np.sqrt(sum_sq / len(time_exp))
   
def normalize_exp_sim(data):
    return np.asarray(data) / max(data)

def shift_time_exp(time_exp, shift):
    time_exp = np.asarray(time_exp) * 1e-3 + shift
    return time_exp


df_ampar = pd.read_csv("../glut_y_compare/collated_errors_AMPAR.csv")

asymmetric_error = np.array(list(zip(df_ampar["lerror"], df_ampar["uerror"]))).T

ax12H.errorbar(df_ampar["y"], df_ampar["AMPAR"], yerr = asymmetric_error, color = 'r', label = "AMPAR Current")
ax12H.set_xlabel("y $\mu m$")
ax12H.set_ylabel("Glutamate & AMPAR")
glut_calc_profile = []
heights = [0.01, 0.5, 1.0, 1.5, 2.0]
for height in heights:
   dendritic_q = glutamate_profile(height, 1)
   if height == heights[0]:
      ref = dendritic_q
   dendritic_q = dendritic_q/ref
   glut_calc_profile.append(dendritic_q)

ax12H.plot(heights, glut_calc_profile, 'b', label = "Glutamate")
ax12H.legend(frameon = False)
glut_rmsd = RMSD_calc(df_ampar["y"], heights, df_ampar["AMPAR"], glut_calc_profile)
ax12H.text(0.25, 0.8, "RMSD: " + str(round(glut_rmsd, 2)), weight = "normal",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax12H.transAxes)

df_exp = pd.read_csv("./Calcium_Helmchen.csv")
norm_exp = normalize_exp_sim(df_exp["conc"])
time_exp = df_exp["time"]
shift_time = 0.5 #simulation pulse comes at 0.5 s
time_exp = shift_time_exp(time_exp, shift_time)
ax21H.plot(np.asarray(time_exp), norm_exp, color = 'r', label = "Exp.")
tag_string = "Source_" + str(source) + "_nb_" + str(num_bursts) + "_ns_" + str(num_spines) + "_ypos_" + str(y_pos) + "_spacing_" + str(spacing) + "_freq_" + str(1.0) + "_CaTau_" + str(CaTau) + "_phi_entire_" + str(phi_entire)
ax21H.set_ylabel("Calcium")
ax21H.set_xlabel("time (s)")

filename = tag_string + "Ca.xml"
tree1 = ET.parse(filename)
tree = [tree1]
sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
norm_Ca = normalize_exp_sim(np.asarray(max_an_y)- min(max_an_y))
shape_dt = 0.001
t_axis = np.linspace(0, len(norm_Ca) * shape_dt, len(norm_Ca))
ax21H.plot( t_axis, norm_Ca, 'b', label = "Sim." ) 
ax21H.legend(frameon = False)
ca_rmsd = RMSD_calc(time_exp, t_axis, norm_exp, norm_Ca)
ax21H.text(0.25, 0.8, "RMSD: " + str(round(ca_rmsd, 2)), weight = "normal",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax21H.transAxes)

df_CaMKII = pd.read_csv("../collated_CaMKII_error.csv")

asymmetric_error = np.array(list(zip(df_CaMKII["lerror"], df_CaMKII["uerror"]))).T

ax22H.errorbar(df_CaMKII["time"], df_CaMKII["life"], color = 'r', yerr = asymmetric_error, label = "Exp")
filename = "Source_32000.0_nb_30_ns_1_ypos_5e-07_spacing_7.5e-07_freq_0.5_CaTau_0.08_phi_entire_6000.0thr286.xml"
tree1 = ET.parse("../CaMKII_quantification/" + filename)
tree = [tree1]
sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
norm_CaMKII = normalize_exp_sim(np.asarray(max_an_y)- min(max_an_y))
shape_dt = 0.2
t_axis = np.linspace(0, len(norm_CaMKII) * shape_dt, len(norm_CaMKII))
t_axis = np.asarray(t_axis) - shift_time
ax22H.plot( t_axis, norm_CaMKII, 'b', label = "Sim." ) 
ax22H.set_xlim([0, 100])
p = patches.Rectangle((0,0), 90, 0.1, linewidth=2, fill=None, hatch='///')
ax22H.add_patch(p)
ax22H.set_xlabel("time (s)")
ax22H.legend(frameon = False)
ax22H.set_ylabel("$CaMKII_{active}$")

camkii_rmsd = RMSD_calc(df_CaMKII["time"], t_axis, df_CaMKII["life"], norm_CaMKII)

ax22H.text(0.25, 0.25, "RMSD: " + str(round(camkii_rmsd, 2)), weight = "normal",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax22H.transAxes)

df_Rac = pd.read_csv("../Rac_quantification/collated_Rac_error.csv")
lerror = np.asarray(df_Rac["lerror"]) * np.asarray(df_Rac["life"])
uerror = np.asarray(df_Rac["uerror"]) * np.asarray(df_Rac["life"])
asymmetric_error = np.array(list(zip(lerror, uerror))).T

ax31H.errorbar(df_Rac["time"], df_Rac["life"], color = 'r', yerr = asymmetric_error, label = "Exp.")
#time_exp = df_Rac["time"]
shift_time = 15 #simulation pulse comes at 15 s
#time_exp = shift_time_exp(np.asarray(time_exp) * 1e3, shift_time)
ax31H.set_ylabel("Rac_GTP")
ax31H.set_xlabel("time (s)")
filename = "Source_32000.0_nb_30_ns_1_ypos_5e-07_spacing_7.5e-07_freq_0.5_CaTau_0.08_phi_entire_6000.0Rac.xml"
tree1 = ET.parse("../Rac_quantification/" + filename)
tree = [tree1]
sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
norm_Rac = normalize_exp_sim(np.asarray(max_an_y)- min(max_an_y))
shape_dt = 0.2
t_axis = np.linspace(0, len(norm_Rac) * shape_dt, len(norm_Rac))
t_axis = np.asarray(t_axis) - shift_time
ax31H.plot( t_axis, norm_Rac, 'b', label = "Sim." ) 
ax31H.set_xlim([0, 250])
ax31H.legend(frameon = False)
p = patches.Rectangle((0,0), 60, 0.1, linewidth=2, fill=None, hatch='///')
ax31H.add_patch(p)

          
rac_rmsd = RMSD_calc(df_Rac["time"], t_axis, df_Rac["life"], norm_Rac)
ax31H.text(0.25, 0.25, "RMSD: " + str(round(rac_rmsd, 2)), weight = "normal",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax31H.transAxes)


ax32H.set_ylabel("IRSp53_act")
ax32H.set_xlabel("time (s)")
filename = "Source_32000.0_nb_30_ns_1_ypos_5e-07_spacing_7.5e-07_freq_0.5_CaTau_0.08_phi_entire_6000.0cytosol.xml"
tree1 = ET.parse("../Rac_quantification/" + filename)
tree = [tree1]
sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
norm_cytosol = normalize_exp_sim(np.asarray(max_an_y)- min(max_an_y))
shape_dt = 0.2
t_axis = np.linspace(0, len(norm_cytosol) * shape_dt, len(norm_cytosol))
t_axis = np.asarray(t_axis) - shift_time
ax32H.plot( t_axis, norm_cytosol, 'b', label = "Sim." ) 
ax32H.set_xlim([0, 250])
ax32H.legend(frameon = False)
p = patches.Rectangle((0,0), 60, 0.1, linewidth=2, fill=None, hatch='///')
ax32H.add_patch(p)
"""
df_diff = pd.read_csv("../DiffC.csv")

for i in range(len(df_diff)):
    print(df_diff["known"][i])
    if df_diff["known"][i] == 1:
        print("Found known molecule")
        color = "r"
    else:
        color = "k"
    ax31H.plot(1/np.cbrt(df_diff["wt"][i]), np.asarray(df_diff["D"][i]), marker = "*", markerfacecolor = color, markeredgecolor = color, linestyle = "None")
ax31H.set_xlabel("1 / $\sqrt[3]{Mol. Wt}$")
ax31H.set_ylabel("D $\mu m^{2}/s$")
ax31H.plot(1 / np.cbrt( np.asarray(df_diff["wt"]) ), np.asarray(df_diff["D"]))
"""
plt.show()


import pandas as pd
import tag
from functools import reduce
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import readXML
import seaborn as sns
import math
import totalEnergy_calc as TE
import wave_speed_calc as Ws
import dendShape
from scipy import optimize



phisat = tag.phisat
mem_diffConst = tag.mem_diffConst
cyt_diffConst = tag.cyt_diffConst
Kb = tag.Kb
part_coeff = TE.part_coeff
print("Kb: ", Kb)

mu0_list = [-2.0, -2.5, -3.0, -3.5, -4.0, -4.5, -5.0, -5.5, -6.0]
kagg_list = [-1.0, -2.0, -4.0, -6.0, -8.0]

def calc_energy(phitot, mu0RT, k_agg):
      mu0 = mu0RT * TE.RT
      conc = (tag.phi_entire - phitot) / (TE.V_cyt * TE.Na)
      ret = optimize.minimize( TE.total_energy, x0=TE.x0, bounds = (TE.rmbounds,TE.tbounds),args=(phitot,TE.phisat, TE.khat, TE.fs, TE.k, 1, 0, 0, 0, conc, mu0, k_agg),method = 'L-BFGS-B', options = {'ftol':1e-12} )
      return ret, conc

def chem_pot_partition(phitot, conc, rm , theta, mu0):
    area_dome = 2 * np.pi * (rm * 1e-6)**2 * (1 - np.cos(theta))
    V_dome = area_dome * TE.thickness
    mem_conc = phitot / (TE.Na * V_dome)
    mu = TE.RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    return (mu + mu0)  #This should give a value in J/molecule.


def mech_potential(hmean, phi, k_agg):
    mismatch_potential = ( 0.5 * TE.khat / phisat ) * (hmean - TE.Cp)**2
    entropy_potential = TE.KBT * (np.log( phi ) - np.log(1 - phi))
    aggregation_potential = TE.KBT * ( 1 + k_agg * phi )
    mechanical_potential = mismatch_potential + entropy_potential + aggregation_potential
    return mechanical_potential

def extract_time_series(hmean_d, theta_d, sstart_d):
    dome_half_width_all = []
    for i in range(len(hmean_d)):
       hmean_list = []
       theta_list = []
       sstart_list = []
       for j in range(1, len(hmean_d.iloc[i])):
           if np.isnan( hmean_d.iloc[i][j] ) == False:
              hmean_list.append( hmean_d.iloc[i][j] )
              theta_list.append( theta_d.iloc[i][j] )
              sstart_list.append( sstart_d.iloc[i][j] )
       x, y, dome_half_width = dendShape.get_shape(sstart_list, theta_list, hmean_list)
       dome_half_width_all.append( dome_half_width )
    return dome_half_width_all 

                   

fig, (ax1, ax2) = plt.subplots(2, 1)
mu_plot_list = []
agg_plot_list = []
speed_list = []
phitot_list = []
kf_list = []
energy_list = []
ph = 1
for mui in mu0_list:
    mu0RT =  mui * (1 / 0.593)
    mu0 = mu0RT * TE.RT
    for k_agg in kagg_list:
        print("mu, kagg: ", mui, k_agg)
        agg_plot_list.append(k_agg)
        mu_plot_list.append(mui)
        tag_string = "phiE_" + str(4000.0) + "_ns_" + str(1) + "_mD_" + str(mem_diffConst) + "_mu0_" + str(mui) + "_spacing_" + str(1e-6) + "_cD_" + str(cyt_diffConst) + "_k_agg_" + str(k_agg) + "_One"
        phi = pd.read_csv(tag_string + "phi.csv")
        hmean = pd.read_csv(tag_string + "hmean.csv")
        theta = pd.read_csv(tag_string + "theta.csv")
        cConc = pd.read_csv(tag_string + "cConc.csv")
        sstart = pd.read_csv(tag_string + "sstart.csv")
        try:
           area = 2 * np.pi * ( 1 / hmean.iloc[-1][ph] )**2 * (1 - np.cos(theta.iloc[-1][ph]))
        except:
           print("File empty - No spines")
           speed_list.append(np.nan)
           phitot_list.append(np.nan)
           kf_list.append(np.nan)
           energy_list.append(np.nan)
        else:   
           Phitot = phi.iloc[-1][ph] * phisat * area
           phitot_list.append(Phitot)
           rm_val = 1 / hmean.iloc[-1][ph]
           theta_val = theta.iloc[-1][ph]
           ret, conc = calc_energy(Phitot, mu0RT, k_agg)
           energy_list.append( ret.fun * 1e-16 / TE.KBT )
           mechanical_potential = mech_potential(hmean.iloc[-1][ph], phi.iloc[-1][ph], k_agg)
           chemical_potential = chem_pot_partition(Phitot, cConc.iloc[-1][ph], rm_val * 1e6, theta_val, mu0)
           print("mu0, kagg: ", mui, k_agg)
           try:
              extract_time_series(hmean, theta, sstart)
           except:
              print("Unstable")
           else:   
              dome_half_width_t = extract_time_series(hmean, theta, sstart)
           print("mu0, k_agg: ", mui, k_agg)
           print("rm, theta, phi: ", rm_val, theta_val, phi.iloc[-1][ph])
           print("cCONC: ", cConc.iloc[-1][ph])
           print("Mech pot, Chem pot: ", mechanical_potential/TE.RT, chemical_potential/TE.RT)
           print("RT: ", TE.RT)
           kf_val = Kb * TE.part_coeff * np.exp( -(mechanical_potential + chemical_potential) / TE.RT )
           kf_list.append( kf_val )
           speed_list.append( Ws.wave_speed_calc( rm_val, theta_val, phi.iloc[-1][ph], cConc.iloc[-1][ph], kf_val )  * 1e9 )

df= pd.DataFrame(columns = ["Agg. coeff.", "mu0", "Speed", "phitot", "Kf", "Energy"])
df["Agg. coeff."] = agg_plot_list
df["mu0"] = mu_plot_list
df["Speed"] = speed_list
df["phitot"] = phitot_list
df["Kf"] = kf_list
df["Energy"] = energy_list

df_phitot = df.pivot("Agg. coeff.", "mu0", "phitot").round(6)
df_speed = df.pivot("Agg. coeff.", "mu0", "Speed").round(6)
df_kf = df.pivot("Agg. coeff.", "mu0", "Kf").round(6)
df_energy = df.pivot("Agg. coeff.", "mu0", "Energy").round(6)

ax2 = sns.heatmap(df_phitot, annot = True, ax = ax2,xticklabels = df_phitot.columns.values.round(2))
ax1 = sns.heatmap(df_energy, annot = True, ax = ax1, xticklabels = df_energy.columns.values.round(2)) 
#ax2 = sns.heatmap(df_kf, annot = True, ax = ax2, xticklabels = df_kf.columns.values.round(2)) 
#ax2.invert_yaxis()
#ax1.invert_yaxis()
ax1.set_title("Energy")
ax2.set_title("$\phi_{tot}$")
plt.tight_layout()        
plt.show()    

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import totalEnergy_wave_calc as TE
import tag
import pandas as pd
import seaborn as sns

search_tolerance = 2e-9
phi_entire = 4000
V_cyt = TE.V_cyt
phisat = TE.phisat
Kb = tag.Kb
diffL = TE.diffL
tot_num_voxels = int(TE.Length * 1e-6 / diffL)
phitot_space = np.linspace(1, phi_entire/2.0, 1000)
mu0_space = [-2.5, -3.5, -4.5,  -5.5, -6.5]
k_agg_list = [-1.0, -2.0, -3.0, -4.0, -6.0, -10.0]


def coeff_calc(phitot, mu0RT, k_agg):
      mu0 = mu0RT * TE.RT
      conc = (phi_entire - phitot) / (TE.V_cyt * TE.Na)
      ret = optimize.minimize( TE.total_energy, x0=TE.x0, bounds = (TE.rmbounds,TE.tbounds),args=(phitot,TE.phisat, TE.khat, TE.fs, TE.k, 1, 0, 0, 0, conc, mu0, k_agg),method = 'L-BFGS-B', options = {'ftol':1e-12} )
      rm_val = ret.x[0] * 1e-6
      theta_val = ret.x[1]
      area = 2 * np.pi * rm_val**2 * (1 - np.cos(theta_val))
      phi = phitot / (area * phisat)
      CI = np.sqrt( 1 / (2 * np.pi * phi * phisat))
      L = rm_val * np.sin( theta_val )
      coeff = 0.5 * phitot**-0.5 * CI * np.sqrt(1 + np.cos(theta_val)) + CI**2 * np.cos(theta_val) / L
      return coeff, rm_val, theta_val, phi, area, conc


def Kf_calc(phitot, rm, phi, conc, mu0RT, theta, k_agg):
    mu0 = mu0RT * TE.RT
    hmean = 1/rm
    chem_potential = TE.chem_pot_partition(phitot, conc, rm * 1e6, theta, mu0)
    mech_potential = TE.mech_potential(hmean, phi, k_agg)
    Kf = Kb * TE.part_coeff * np.exp( -(chem_potential + mech_potential) / TE.RT )
    return Kf

def wave_speed(Kf, coeff, phitot, rm, theta):
    num_cyt_per_voxel = (phi_entire - phitot) / tot_num_voxels
    L = rm * np.sin(theta)
    num_voxels_under_dome = 2 * L / diffL
    num_mem_per_voxel = phitot / num_voxels_under_dome
    return coeff * (Kf * num_cyt_per_voxel - Kb * num_mem_per_voxel)

def crossing_phitot(mu0RT, k_agg):
    ws = []
    for pt in phitot_space:
       coeff, rm, theta, phi, area, conc = coeff_calc(pt, mu0RT, k_agg)
       Kf = Kf_calc(pt, rm, phi, conc, mu0RT, theta, k_agg)
       ws.append( wave_speed(Kf, coeff, pt, rm, theta) )
    if k_agg == select_kagg:   
      plt.figure(101)
      plt.plot(phitot_space, ws, label = "mu0 = " + str(mu0RT * 0.593))
      plt.legend()
      plt.xlabel("$\phi_{tot}$")
      plt.ylabel("Wave speed m/s")

    smallest = 1e-6
    for w in range(len(ws)):
        if np.abs(ws[w]) < smallest:
            smallest = ws[w]
            index = w
            
    near_zero_speed = smallest
    cross_phitot = phitot_space[index]
    print("WS, phitot: ", near_zero_speed, cross_phitot)
    return near_zero_speed, cross_phitot

select_kagg = -10
mu_list = []
agg_list = []
speed_list = []
stable_phitots = []
for mui in mu0_space:
  mu0RT = mui * (1 / 0.593)
  for agg in k_agg_list:
      mu_list.append(mui)
      agg_list.append(agg)
      stable_speed, stable_phitot = crossing_phitot(mu0RT, agg)
      if np.abs(stable_speed) > search_tolerance:
          speed_list.append( stable_speed * 1e9 )
          stable_phitots.append( 0 )
      else:    
          speed_list.append( stable_speed * 1e9 )
          stable_phitots.append( stable_phitot )
df = pd.DataFrame( columns = ["Agg. coeff.", "mu0", "Speed", "Phitot"])
df["Agg. coeff."] = agg_list
df["mu0"] = mu_list
df["Speed"] = speed_list
df["Phitot"] = stable_phitots

fig, (ax1, ax2) = plt.subplots(2,1)
df_speed = df.pivot("Agg. coeff.", "mu0", "Speed").round(6)
df_phitot = df.pivot("Agg. coeff.", "mu0", "Phitot").round(6)
ax1 = sns.heatmap(df_speed, annot = True, ax = ax1, xticklabels = df_speed.columns.values.round(2))
ax2 = sns.heatmap(df_phitot, annot = True, ax = ax2, xticklabels = df_phitot.columns.values.round(2))
ax1.set_title("Wave speed (nm/s)")
ax2.set_title("$\phi_{tot}$")

plt.tight_layout()
print(df)
plt.show()      


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_context("poster")


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.text(-0.1, 1.1, 'A', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax1.transAxes)
ax2.text(-0.1, 1.1, 'B', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax2.transAxes)
ax3.text(-0.1, 1.1, 'C', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax3.transAxes)
ax4.text(-0.1, 1.1, 'D', weight = "bold",
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax4.transAxes)

Length = 20e-6
dendR = 0.5e-6
Na = 6.0223e23
V_cyt = np.pi * dendR**2 * Length

k_agg_list = [-45.0, -50.0, -55.0, -60.0, -65.0]

phi_entire_list = [2000.0, 4000.0, 6000.0]

sample_k_agg = -45.0
sample_phi_entire = 2000.0

def plot_error(phi_entire, axis):
    dE = []
    directory = "../Fig2_new/"
    k_agg_list = [-45.0, -50.0, -55.0, -60.0]
    for k in k_agg_list:
       df_wave = pd.read_csv(directory + "/wave_sol/minData" + "_" + str(k) + ".csv")
       wave_energy = list(df_wave[df_wave["phi_entire"] == phi_entire]["energy"])
       df_sim = pd.read_csv(directory + "/stability_comparison_L20/sim.csv")
       df_sim_energy = df_sim[df_sim["phi_entire"] == phi_entire]
       df_sim_energy = df_sim_energy[df_sim_energy["Aggregation"] == k]
       sim_energy = list(df_sim_energy["energy"])
       dE.append(abs( (wave_energy[0] - sim_energy[0]) / wave_energy[0] ) * 100)
       print("Wave energy")
       print(wave_energy)
       print("Sim energy")
       print(sim_energy)
    conc = 1e3 * phi_entire / (Na * V_cyt)
    axis.plot(k_agg_list, dE, marker = "o", label = "Conc: " + str(round(conc, 2)) + " $\mu M$")
    axis.set_ylabel("Error")
    axis.set_xlabel("$k_{agg}$")
    axis.set_xticks(k_agg_list)
    axis.invert_xaxis()
    axis.legend(frameon = False)
    axis.spines[['right', 'top']].set_visible(False)

def k_agg_dE(phi_entire, axis):
    dE = []
    for k in k_agg_list:
        file_name = "mu-4.5_probability_" + str(k) + "_" + str(phi_entire) + "_5.0nm.csv"
        df = pd.read_csv(file_name)
        df_highEnergy = df["highEnergy"]
        min_highEnergy = np.nanmin(df_highEnergy)
        df_lowEnergy = df["lowEnergy"]
        min_lowEnergy = np.nanmin(df_lowEnergy)
        dE.append(min_highEnergy - min_lowEnergy)
    conc = 1e3 * phi_entire / (Na * V_cyt)
    axis.plot(k_agg_list, np.asarray(dE) / 1e3, label = str(round(conc, 2)) + " $\mu M$" )
    axis.invert_xaxis()
    axis.set_ylabel("dE ($10^{3} K_{B}T$)")
    axis.set_xlabel("$k_{agg}$")
    axis.legend(frameon = False)
    axis.spines[['right', 'top']].set_visible(False)

def sample_plot(axis):
        file_name = "mu-4.5_probability_" + str(sample_k_agg) + "_" + str(sample_phi_entire) + "_5.0nm.csv"
        df = pd.read_csv(file_name)
        df_highEnergy = df["highEnergy"]
        phitot_highE = df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["phitot"]
        df_lowEnergy = df["lowEnergy"]
        phitot_lowE = df[df["lowEnergy"] == np.nanmin(df["lowEnergy"])]["phitot"]
        axis.plot(df["phitot"], np.asarray(df_highEnergy) / 1e3, "k")
        axis.plot(df["phitot"], np.asarray(df_lowEnergy) / 1e3, color = "k",  linestyle = "-.")
        axis.set_ylabel("E ($10^{3} K_{B}T$)")
        axis.set_ylim(np.nanmin(df_highEnergy) / 1e3, 0)
        axis.plot(phitot_lowE, np.nanmin(df_lowEnergy) / 1e3, marker = "o", mfc = "b", mec = "b")
        axis.plot(phitot_highE, np.nanmin(df_highEnergy) / 1e3, marker = "o", mfc = "b", mec = "b")
        axis.text(phitot_highE, np.nanmin(df_highEnergy) / 1e3 + 0.1, "$E_{Sharp}^{min}$")
        axis.text(phitot_lowE, np.nanmin(df_lowEnergy) / 1e3 + 0.1, "$E_{Shallow}^{min}$")
        axis.legend(frameon = False)
        axis.set_xlabel("$\phi_{tot}$")
        axis.spines[['right', 'top']].set_visible(False)

   
sample_plot(ax2)   


for phi_entire in phi_entire_list:
    k_agg_dE(phi_entire, ax3)

for phi_entire in phi_entire_list:
    plot_error(phi_entire, ax4)

ax1.spines[['right', 'top', 'bottom', 'left']].set_visible(False)
ax1.set_xticks([])
ax1.set_yticks([])

plt.show()


import numpy as np
import pandas as pd
import kf_calc

phi_entire = 6000
Length = 20
dendDia = 1e-6
thickness = 5e-9
V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6
Kb = 40.0
Na = 6.022140e23
k_agg = -45.0
mu0_exp = -4.5
read = input("File name")
df = pd.read_csv(read)
he = list(df["highEnergy"])
pe = list(df["phitot"])
print("Min: ", np.nanmin(he))
min_index = he.index(np.nanmin(he))
print("Min phitot: ", pe[min_index])
print(df.iloc[he.index(min(he))]["highEnergy"])
rm = df.iloc[he.index(min(he))]["minRm"]
theta = df.iloc[he.index(min(he))]["minTheta"]
rp = df.iloc[he.index(min(he))]["minRp"]
phitot = df.iloc[he.index(min(he))]["phitot"]
phi = df.iloc[he.index(min(he))]["maxHPHI"]

print("Min: ", rm, theta, rp)
print("PHI: ", phi)

rm = df.iloc[he.index(min(he))]["maxRm"]
theta = df.iloc[he.index(min(he))]["maxTheta"]
rp = df.iloc[he.index(min(he))]["maxRp"]
print("Max: ", rm, theta, rp)

Kf, chem_potential, mech_potential = kf_calc.Kf_calc(phi_entire, phitot, rm, phi, mu0_exp, theta, k_agg)

def flux_num(Kf, phitot, rm, theta):
    L = 2 * rm * np.sin(theta)
    num_div = Length / L
    volume_under_dome = V_cyt / num_div
    num_cyt_molecules_under_dome = (phi_entire - phitot) / num_div
    return ( 1 / (volume_under_dome * Na) ) * (Kf * num_cyt_molecules_under_dome - Kb * phitot)

def conc_dome_calc(phitot, rm, theta):
    area = 2 * np.pi * (rm)**2 * (1 - np.cos(theta))
    V_dome = area * thickness
    mem_conc = phitot / (Na * V_dome)
    return mem_conc

def conc_dend_calc(phitot, phi_entire):    
    conc = (phi_entire - phitot) / ( V_cyt * Na )
    return conc

def flux_calc(dome, dend, Kf, Kb):
    return Kf * dend - Kb * dome

print("Kf: ", Kf)
conc_dome = conc_dome_calc(phitot, rm, theta)
conc_dend = conc_dend_calc(phitot, phi_entire)
print(conc_dome, conc_dend)
flux = flux_calc(conc_dome, conc_dend, Kf, Kb)
print("Flux: ", flux)
fluxN = flux_num(Kf, phitot, rm, theta)
print("Flux Num: ", fluxN)

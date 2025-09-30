from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint , NonlinearConstraint
from multiprocessing import Pool
from functools import partial

import numpy as np
import pandas as pd
import random
from scipy.optimize import Bounds

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--k_agg1")
parser.add_argument("--phi_entire")
args = parser.parse_args()

KBT = 1.380649e-23 * 300
Na = 6.022140e23
RT = 2479 / Na
RT_kcal_mol = 0.593 #kcal/mol
k = 60 * KBT
khat = 40 * KBT
kentropy = KBT
thickness = 5e-9
fm = 0.0e-6
dendDia = 1e-6
Rc = 0.5 * dendDia
k_bind_neq = 0.0
k_bind_eq = 0.0
khat_KBT = khat/KBT
phisat =  1/50e-18
fs = 55 * 1e-7
Cp = 1 / 20e-9
f_cm = 1
Length = 20
xtot = (Length) * 1e-6       
diffL = 20e-9
zero_sub = 1.0
#x0 = [0.55, 0.45, 0.35]  #Note that the init. guess of Rp should be <=0.5 for consistency across all initial guesses.
#x0 = [0.5, 0.5, 0.5]  #Note that the init. guess of Rp should be <=0.5 for consistency across all initial guesses.
albounds = (0.001,1.5)
rpbounds = (0.001,5.0)
tbounds = (0.01, 1.5)
bounds = Bounds([0.001, 0.1, 0.001], [5.0, 1.5, 5.0]) 
nbounds = (1.0, 10.0)
part_coeff = 1
mu0_exp = -4.5
print(args.k_agg1)
k_agg = float(args.k_agg1)
k_agg2 = -k_agg/2.0
k_agg3 = 0.0
k_agg4 = 0.0
phi_entire = float(args.phi_entire)
V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6
#mu0_exp = -8.6 #kcal/mol Takemura Nat. Sci. Rep., 2017
print("RT: ", RT)
print("KBT: ", KBT)

print("Arguments: ", k_agg, k_agg2)

def calc_criterion(conc):
    const1 = (thickness * Na * conc) / phisat
    mu0 = calc_mu0()
    const2 = np.abs(mu0) / RT
    return const1 * np.exp(const2)


def mech_potential(hmean, phi):
    mismatch_potential = ( 0.5 * khat / phisat ) * (hmean - Cp)**2
    entropy_potential = KBT * ( np.log( phi ) - np.log( 1 - phi ) )
    agg_potential = 1 + k_agg * phi + k_agg2 * 3 * phi**2
    aggregation_potential = KBT * agg_potential
    mechanical_potential = mismatch_potential + entropy_potential + aggregation_potential
    return mechanical_potential    

def chem_potential(phitot, conc, area, mu0):
    V_dome = area * thickness
    mem_conc = phitot / (Na * V_dome)
    mu = RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    return mu + mu0

def calc_sagit(alength, theta, rp):
    rm = alength / theta
    r0 = (rm + rp) * np.sin(theta)
    eq_coeff = [1, -2 * Rc, 0.25 * (2 * r0)**2]
    sol = np.roots(eq_coeff)
    sagitta = np.real(sol)[1]
    return sagitta

def calc_dome_area(alength, theta):
    return 2 * np.pi * (alength / theta)**2 * (1 - np.cos(theta))

def calc_saddle_area(alength, theta, rp):
    rm = alength / theta
    r0 = (rp + rm) * np.sin(theta)
    return 2 * np.pi * rp * ( r0 * theta + rp * ( np.cos(theta) - 1 ) )

def calc_mu0():
    mu0RT = mu0_exp / RT_kcal_mol
    mu0 = mu0RT * RT
    return mu0

def calc_phi(xp): 
    area = calc_dome_area(xp[0] * 1e-6, xp[1])
    return phitot / (area * phisat)

def chem_pot_partition(phitot, conc, alength , theta, mu0, area):
    rm = alength / theta
    V_dome = area * thickness
    mem_conc = phitot / (Na * V_dome)
    #print("Area: ", area)
    #print("MEM CONC: ", mem_conc)
    mu = RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    phi = phitot / (area * phisat)
    return (mu + mu0) * phitot  #This should give a value in J/molecule.


def aggregation_energy(ph_en, phitot, area):
    polynomial_eval = ph_en + 0.5 * k_agg * ph_en**2 + k_agg2 * ph_en**3 + k_agg3 * ph_en**4 + k_agg4 * ph_en**5
    #polynomial_eval = 0.5 * k_agg * ph_en**2 
    return area * (kentropy * phisat) * polynomial_eval


def entropy( ph_en, phitot, area, khat):
    entropy = area * (kentropy * phisat) * (ph_en * np.log(ph_en) + (1 - ph_en) * np.log(1 - ph_en))
    return entropy

def mismatch( alength, theta, phi, phitot, area ):
    # Note that the  area term cancels out with the phitot/(area)
    #return area * 0.5 * khat * phi * (2/rm - Cp) * (2/rm - Cp)
    rm = alength / theta
    return 0.5 * khat * (phitot / phisat) * (1/rm - Cp) * (1/rm - Cp)

def fm_tension( ph_en, phitot, area ):
    return area * fs

def fp_tension(alength, theta, rp):
    rm = alength / theta
    saddle_area = calc_saddle_area(alength, theta, rp)
    return saddle_area * fs

def Fplus(alength, theta, rp):
    rm = alength / theta
    r0 = (rp + rm) * np.sin(theta)
    r_psi = rm * np.sin(theta)
    saddle_area = calc_saddle_area(alength, theta, rp)
 

    mean_bend_energy = saddle_area * 0.5 * k * ( ( 1/rp + 1/r_psi ) / 2.0 )**2
    gauss_bend_energy = saddle_area * 0.5 * k * 1 / (rp * r_psi)

    mean_energy = fs * np.pi * r0**2
    #total_bend_energy = mean_bend_energy + gauss_bend_energy - mean_energy
    total_bend_energy = mean_bend_energy - mean_energy + gauss_bend_energy
    return total_bend_energy

def Fzero(alength, theta, rp):
    rm = alength / theta
    r0 = (rm + rp) * np.sin(theta)
    return fs * np.pi * r0**2

def constraint_abs(alength, theta):
    rm = alength / theta
    return fs * np.pi * rm**2 * (1 - np.cos(theta))**2

#def constraint_abs(rm, theta):  #Included twice the tension term. So no need to add the tension energy in the total_energy
#    return fs * np.pi * rm**2 * (1 - np.cos(theta)) * (3 - np.cos(theta))

def constraint_ratio(alength, theta):
    rm = alength / theta
    return fs * np.pi * rm**2 * (1 - np.cos(theta))**3 / (1 + np.cos(theta))

def dome_bending(alength, theta, ph_en, area):
    rm = alength / theta
    return area * (k / (2 * rm * rm))

def calc_dome_volume(alength, theta):
    rm = alength / theta
    sphere_rad = rm * np.sin(theta)
    return (2 * np.pi / 3.0) * sphere_rad**3

def calc_saddle_volume(alength, theta, rp):
    rm = alength / theta
    r0 = (rm + rp) * np.sin(theta)
    first_term = r0**2 * rp * (1 - np.cos(theta))
    second_term = r0 * rp**2 * (theta - np.sin(theta) * np.cos(theta))
    third_term = (rp**3 / 12.0) * (np.cos(3 * theta) - 9 * np.cos(theta) + 8)
    volume = np.pi * (first_term + second_term + third_term)
    return volume

def pv_work(alength, theta, rp):
    rm = alength / theta
    cyl_rad = 0.5 * dendDia
    r0 = (rm + rp) * np.sin(theta)
    #cyl_volume = np.pi * cyl_rad**2 * rm * np.sin(theta)
    cyl_volume = np.pi * cyl_rad**2 * r0
    dome_volume = calc_dome_volume(rm, theta)
    saddle_volume = calc_saddle_volume(rm, theta, rp)
    total_volume = dome_volume + saddle_volume
    volume = -cyl_volume + total_volume
    #volume = total_volume
    pressure = 2 * fs / cyl_rad
    #pressure = 2 * fs / r0
    return pressure * volume

def constraint_sagitt(alength , theta, rp):
    rm = alength / theta
    saddle_area = calc_saddle_area(alength, theta, rp)
    dome_area = calc_dome_area(alength, theta)
    total_area = saddle_area + dome_area
    deltaH = calc_sagit(alength, theta, rp)
    return (deltaH / dendDia) * fs * total_area

def constraint_yR(alength, theta, rp):
    rm = alength / theta
    r0 = (rm + rp) * np.sin(theta)
    saddle_area = calc_saddle_area(alength, theta, rp)
    dome_area = calc_dome_area(alength, theta)
    total_area = saddle_area + dome_area
    gamma = np.arctan(r0 / Rc)
    yR = Rc / np.cos(gamma)
    cost_var = (yR - Rc)**2 / Rc**2
    #cost_var = (yR - Rc)**3 / Rc**3
    return cost_var * total_area * fs


def area_constraint(alength, theta, rp):
    rm = alength / theta
    r0 = (rp + rm) * np.sin( theta )
    disk_area = np.pi * r0 * r0
    return fs * disk_area * (2 * rm / dendDia )

def total_energy(x, phitot):
    alength = x[0] * 1e-6
    theta = x[1]
    rp = x[2] * 1e-6
    area = calc_dome_area(alength, theta)
    ph_en = phitot / (area * phisat)
    mismatch_E = mismatch( alength, theta, ph_en, phitot, area )
    entropy_E = entropy( ph_en, phitot, area, khat)
    bending_E = dome_bending(alength, theta, ph_en, area)
    tension_E = fm_tension( ph_en, phitot, area )
    agg_en = aggregation_energy(ph_en, phitot, area)
    mu0 = calc_mu0()
    conc = (phi_entire - phitot) / ( V_cyt * Na )
    chem_pot = chem_pot_partition(phitot, conc, alength ,  theta, mu0, area)
    
    saddle_tension = fp_tension(alength, theta, rp)
    saddle_bend = Fplus(alength, theta, rp)
    zero_surface = Fzero(alength, theta, rp)
    cross_section = np.pi * (0.5 * dendDia)**2
    pressure_volume_work = pv_work(alength, theta, rp)
    #constraint_energy = constraint_sagitt(alength, theta, rp)
    constraint_rad_gap =  constraint_yR(alength, theta, rp)
    #print(constraint_rad_gap)
    constraint_area_energy = area_constraint(alength, theta, rp)
    #total_energy = mismatch_E + entropy_E + bending_E + tension_E + agg_en + chem_pot + saddle_tension + saddle_bend + pressure_volume_work
    total_energy = mismatch_E + entropy_E + bending_E + tension_E + agg_en + chem_pot + saddle_tension + saddle_bend + constraint_rad_gap
    return total_energy / KBT


linear_constraint1 = LinearConstraint( [[1, 0, 0]], [0.001], [5.0] )
linear_constraint2 = LinearConstraint( [[0, 1, 0]], [0.001], [1.5] )
linear_constraint3 = LinearConstraint( [[0, 0, 1]], [0.001], [5.0] )
phitot_space = np.linspace(100, phi_entire/1.0, 25)
df_mins = pd.read_csv("./mu" + str(mu0_exp) + "_probability_" + str(k_agg) + "_" + str(phi_entire) + "_5.0nm.csv")
phitot_min_lowEnergy = int(df_mins[df_mins["lowEnergy"] == np.nanmin(df_mins["lowEnergy"])]["phitot"])
phitot_min_highEnergy = int(df_mins[df_mins["highEnergy"] == np.nanmin(df_mins["highEnergy"])]["phitot"])
phitot_space = list(phitot_space) + [phitot_min_lowEnergy] + [phitot_min_highEnergy]

def gen_random():
    total_samples = []
    for n in range(n_opt_it):
        dx0 = []
        for i in range(3):
              dx0.append(random.uniform(0.001, 1.0))
        total_samples.append(np.asarray(dx0))
    return total_samples

def split_samples(init_samples, num_split):
    samples = []
    start = 0
    each_set_n = int(len(init_samples) / num_split)
    for nsplit in range(num_split):
       sample_set = [] 
       for isamp in range(start, start + each_set_n):
           sample_set.append(init_samples[isamp])
           start = isamp
       samples.append( sample_set )
    return samples   

def calc_params(opt_data, phitot):
   print(opt_data) 
   theta = opt_data.x[1]
   alength = opt_data.x[0] * 1e-6
   rm = alength / theta
   rp = opt_data.x[2] * 1e-6
   r0 = (rm + rp) * np.sin(theta)
   area = calc_dome_area(alength, theta)
   ph_it = phitot / (area * phisat)
   return rm, theta, rp, area, ph_it, alength, r0

def calc_energy(init_con, phitot, conc):
   phi_Llim = 0.0
   #phi_Hlim = calc_criterion(conc) + epsilon
   phi_Hlim = 1.0
   Limits = [phi_Llim, phi_Hlim]

   def calc_phi(xp):
       area = calc_dome_area( (xp[0]/xp[1]) * 1e-6, xp[1])
       phi = phitot / (area * phisat)
       return phi
   nonlinear_constraint = NonlinearConstraint(calc_phi, Limits[0], Limits[1])
   try:
      ret = minimize(total_energy, init_con, method='trust-constr', args = (phitot),  constraints = [linear_constraint1, linear_constraint2, linear_constraint3, nonlinear_constraint])
   except:
      print("Contains NaN") 
      ret = []
   return ret

def calc_failSafe():
    init_samples = gen_random()
    return init_samples

x0 = np.array([0.5, 0.5, 0.5]) 
n_opt_it = 100

def find_branch_probability(phitot, conc):
    samples = split_samples(init_samples, 4)
    lowCount = 0
    highCount = 0
    for samp in samples:
       with Pool() as pool:
           result = pool.map(partial(calc_energy, phitot = phitot, conc = conc), samp)
           for res in result:
              if res != []:
                rm, theta, rp, area, ph_it, alength, r0 = calc_params(res, phitot)
                if ph_it < 0.1 and lowCount == 0:
                     low_rm = rm
                     low_theta = theta
                     low_rp = rp
                     lowCount = lowCount + 1
                if ph_it > 0.1 and highCount == 0:
                     print("High phi found")
                     high_rm = rm
                     high_theta = theta
                     high_rp = rp
                     highCount = highCount + 1
    try:
       print(low_rm, low_theta, low_rp)
       state = "1"
    except:
       print("Lower branch not found")
       low_rm = 1
       low_theta = 1
       low_rp = 1
       state = "0"
    try:
       print(high_rm, high_theta, high_rp)
       state = state + "1"
    except:
       print("Higher branch not found")
       high_rm = 1
       high_theta = 1
       high_rp = 1
       state = state + "0"
    return low_rm, low_theta, low_rp, high_rm, high_theta, high_rp, state

init_samples = calc_failSafe()    
lowRm_list = []
lowTheta_list = []
lowRp_list = []
lowME = []
lowAgg = []
lowEn = []
lowBE = []
lowTE = []
lowCE = []
lowConE = []
lowMechPot = []
lowTotal = []

highME = []
highAgg = []
highEn = []
highBE = []
highTE = []
highCE = []
highConE = []
highMechPot = []
highTotal = []

highRm_list = []
highTheta_list = []
highRp_list = []

state_list = []

def individual_energies(x, phitot):
    alength = x[0] * 1e-6
    theta = x[1]
    rp = x[2] * 1e-6
    rm = alength / theta
    mu0 = calc_mu0()
    conc = (phi_entire - phitot) / ( V_cyt * Na )
    area = calc_dome_area(alength, theta)
    ph_en = phitot / (area * phisat)
    mismatch_E = mismatch( alength, theta, ph_en, phitot, area ) / KBT
    entropy_E = entropy( ph_en, phitot, area, khat) / KBT
    bending_E = dome_bending(alength, theta, ph_en, area) / KBT
    tension_E = fm_tension( ph_en, phitot, area ) / KBT
    agg_en = aggregation_energy(ph_en, phitot, area) / KBT
    chem_pot = chem_pot_partition(phitot, conc, alength ,  theta, mu0, area) / KBT
    mech_pot = mech_potential(1/rm, ph_en) / KBT
    constraint_rad_gap =  constraint_yR(alength, theta, rp) / KBT
    return mismatch_E, entropy_E, bending_E, tension_E, agg_en, chem_pot, constraint_rad_gap, mech_pot


for phitot in phitot_space:
   conc = (phi_entire - phitot) / ( V_cyt * Na )
   low_rm, low_theta, low_rp, high_rm, high_theta, high_rp, state = find_branch_probability(phitot, conc)
   lowRm_list.append(low_rm)
   lowTheta_list.append(low_theta)
   lowRp_list.append(low_rp)
   mismatch_E, entropy_E, bending_E, tension_E, agg_en, chem_pot, constraint_rad_gap, mech_pot = individual_energies([low_rm * 1e6 * low_theta, low_theta, low_rp * 1e6],  phitot)
   Etotal = total_energy([low_rm * 1e6 * low_theta, low_theta, low_rp * 1e6], phitot)
   lowME.append(mismatch_E)
   lowAgg.append(agg_en)
   lowEn.append(entropy_E)
   lowBE.append(bending_E)
   lowTE.append(tension_E)
   lowCE.append(chem_pot)
   lowMechPot.append(mech_pot)
   lowConE.append(constraint_rad_gap)
   lowTotal.append(Etotal)
   
   highRm_list.append(high_rm)
   highTheta_list.append(high_theta)
   highRp_list.append(high_rp)
   state_list.append(state)
   mismatch_E, entropy_E, bending_E, tension_E, agg_en, chem_pot, constraint_rad_gap, mech_pot = individual_energies([high_rm * 1e6 * high_theta, high_theta, high_rp * 1e6],  phitot)
   Etotal = total_energy([high_rm * 1e6 * high_theta, high_theta, high_rp * 1e6], phitot)
   highME.append(mismatch_E)
   highAgg.append(agg_en)
   highEn.append(entropy_E)
   highBE.append(bending_E)
   highTE.append(tension_E)
   highCE.append(chem_pot)
   highConE.append(constraint_rad_gap)
   highMechPot.append(mech_pot)
   highTotal.append(Etotal)
   print("Next phitot: ", phitot)

df = pd.DataFrame()
df["phitot"] = phitot_space
df["LRm"] = lowRm_list
df["LTheta"] = lowTheta_list
df["LRp"] = lowRp_list

df["LME"] = lowME
df["LBE"] = lowBE
df["LEN"] = lowEn
df["LTE"] = lowTE
df["LCE"] = lowCE
df["LMPot"] = lowMechPot
df["LCON"] = lowConE
df["LAG"] = lowAgg
df["lowTotal"] = lowTotal

df["HRm"] = highRm_list
df["HTheta"] = highTheta_list
df["HRp"] = highRp_list

df["HME"] = highME
df["HBE"] = highBE
df["HEN"] = highEn
df["HTE"] = highTE
df["HCE"] = highCE
df["HMPot"] = highMechPot
df["HCON"] = highConE
df["HAG"] = highAgg
df["highTotal"] = highTotal
df["state"] = state_list

df.to_csv("./" + str(args.phi_entire) + "_" +str(args.k_agg1) + "geom_param.csv")



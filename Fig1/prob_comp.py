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
diffL = 2e-9
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

def mech_potential(alength, theta, phitot):
    rm = alength / theta
    hmean = 1/rm
    area = calc_dome_area(alength, theta)
    phi = phitot / (phisat * area)
    mismatch_potential = ( 0.5 * khat / phisat ) * (hmean - Cp)**2
    entropy_potential = KBT * ( np.log( phi ) - np.log( 1 - phi ) )
    agg_potential = 1 + k_agg * phi + k_agg2 * 3 * phi**2
    aggregation_potential = KBT * agg_potential
    print("Area: ", area)
    print("Phisat: ", phisat)
    print("PHI: ", phi)
    print("Ind. poten: ", mismatch_potential, entropy_potential, aggregation_potential)
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
    mu =  RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    #xt = phi_entire / Na
    #xm = phitot / Na
    #xc = (phi_entire - phitot) / Na 
    #mu =  RT * np.log( (xm + 1e-24) / (xc + 1e-24) )
    #phi = phitot / (area * phisat)
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
    mech_pot = mech_potential(alength, theta, phitot) * phitot
    #print("MECH POT: ", mech_pot)
    
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
phitot_space = np.linspace(100, phi_entire, 75)

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
n_opt_it = 50

def find_branch_probability(phitot, conc):
    low_phi_list = []
    high_phi_list = []
    low_energy_list = []
    high_energy_list = []
    
    high_CP_list = []
    high_MP_list = []
    """
    min_low_en = []
    max_low_en = []
    min_high_en = []
    max_high_en = []
    min_high_CP = []
    max_high_CP = []
    min_high_MP = []
    max_high_MP = []
    """
    samples = split_samples(init_samples, 10)
    for samp in samples:
       with Pool() as pool:
           result = pool.map(partial(calc_energy, phitot = phitot, conc = conc), samp)
           for res in result:
              if res != []:
                rm, theta, rp, area, ph_it, alength, r0 = calc_params(res, phitot)
                mu0 = calc_mu0()
                chem_pot_temp = chem_pot_partition(phitot, conc, rm * theta , theta, mu0, area)
                mech_pot_temp =  mech_potential(rm * theta, theta, phitot) * phitot
                print("Mech pot now: ", mech_pot_temp, rm, theta, phitot, ph_it)
                if ph_it < 0.1:
                     low_phi_list.append(ph_it)
                     low_energy_list.append(res.fun)
                if ph_it > 0.1:
                     high_phi_list.append(ph_it)
                     high_energy_list.append(res.fun)
                     high_CP_list.append(chem_pot_temp)
                     high_MP_list.append(mech_pot_temp)
                     print("PHI from calc: ", ph_it)
                     print("high MP list: ", high_MP_list)


    total_success = len(low_phi_list) + len(high_phi_list)       
    if total_success != 0.0:
       low_phi_p = len(low_phi_list) / total_success
       high_phi_p = len(high_phi_list) / total_success
       try:
         low_energy_av = sum(low_energy_list) / len(low_energy_list)
         #low_energy_av = low_energy_list[0]
         min_low_en = min(low_energy_list)
         max_low_en = max(low_energy_list)
       except:
         low_energy_av = np.nan
         min_low_en = np.nan
         max_low_en = np.nan
       try:  
         high_energy_av = sum(high_energy_list) / len(high_energy_list)
         #high_energy_av = high_energy_list[0]
         min_high_en = min(high_energy_list)
         max_high_en = max(high_energy_list)
         min_high_CP = min(high_CP_list)
         max_high_CP = max(high_CP_list)
         min_high_MP = min(high_MP_list)
         max_high_MP = max(high_MP_list)
       except:
         high_energy_av = np.nan
         min_high_en = np.nan
         max_high_en = np.nan
         min_high_CP = np.nan
         max_high_CP = np.nan
         min_high_MP = np.nan
         max_high_MP = np.nan
    else:   
       print("None passed") 
       low_phi_p = np.nan
       high_phi_p = np.nan
       low_energy_av = np.nan
       high_energy_av = np.nan
    try:
        maxLPHI = max(low_phi_list)
    except:
        maxLPHI = 0
    try:
        maxHPHI = max(high_phi_list)
    except:
        maxHPHI = 0
    return low_phi_p, high_phi_p, low_energy_av, high_energy_av, min_low_en, max_low_en, min_high_en, max_high_en, maxLPHI, maxHPHI, min_high_CP, max_high_CP, min_high_MP, max_high_MP

init_samples = calc_failSafe()    
low_phi_p_list = []
high_phi_p_list = []
low_energy_list = []
high_energy_list = []
minLE_list = []
maxLE_list = []
minHE_list = []
maxHE_list = []
maxLPHI_list = []
maxHPHI_list = []
minCP_list = []
maxCP_list = []
minMP_list = []
maxMP_list = []
for phitot in phitot_space:
   conc = (phi_entire - phitot) / ( V_cyt * Na )
   low_phi_p, high_phi_p, lowE, highE, minLE, maxLE, minHE, maxHE, maxLPHI, maxHPHI, min_high_CP, max_high_CP, min_high_MP, max_high_MP = find_branch_probability(phitot, conc)
   low_phi_p_list.append(low_phi_p)
   high_phi_p_list.append(high_phi_p)
   low_energy_list.append(lowE)
   high_energy_list.append(highE)
   minLE_list.append(minLE)
   maxLE_list.append(maxLE)
   minHE_list.append(minHE)
   maxHE_list.append(maxHE)
   maxLPHI_list.append(maxLPHI)
   maxHPHI_list.append(maxHPHI)
   minCP_list.append(min_high_CP)
   maxCP_list.append(max_high_CP)
   minMP_list.append(min_high_MP)
   maxMP_list.append(max_high_MP)
   print("Next phitot: ", phitot)
   print("Mech pot: ", min_high_MP)
   print("Mech pot: ", max_high_MP)
   print(low_phi_p_list)
df = pd.DataFrame()
df["probLow"] = low_phi_p_list
df["probHigh"] = high_phi_p_list
df["lowEnergy"] = low_energy_list
df["minLE"] = minLE_list
df["maxLE"] = maxLE_list
df["minHE"] = minHE_list
df["maxHE"] = maxHE_list
df["maxLPHI"] = maxLPHI_list
df["maxHPHI"] = maxHPHI_list
df["highEnergy"] = high_energy_list
df["minCP"] = minCP_list
df["maxCP"] = maxCP_list
df["minMP"] = minMP_list
df["maxMP"] = maxMP_list
df["phitot"] = phitot_space

df.to_csv("./mu" + str(mu0_exp) + "_probability_" + args.k_agg1 + "_" + args.phi_entire + "_" + str(thickness * 1e9) + "nm.csv")


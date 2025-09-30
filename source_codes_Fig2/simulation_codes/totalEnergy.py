from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from multiprocessing import Pool
from functools import partial

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
import pandas as pd
import random
from scipy.optimize import Bounds
sns.set_context("poster")

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
x0 = [0.5, 0.5, 0.5]  #Note that the init. guess of Rp should be <=0.5 for consistency across all initial guesses.
albounds = (0.001,1.5)
rpbounds = (0.001,5.0)
tbounds = (0.01, 1.5)
bounds = Bounds([0.001, 0.1, 0.001], [5.0, 1.5, 5.0]) 
nbounds = (1.0, 10.0)
part_coeff = 1
phi_entire = 4000
V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6
print("RT: ", RT)
print("KBT: ", KBT)

def mech_potential(hmean, phi, k_agg, k_agg2):
    mismatch_potential = ( 0.5 * khat / phisat ) * (hmean - Cp)**2
    entropy_potential = KBT * ( np.log( phi ) - np.log( 1 - phi ) )
    agg_potential = 1 + k_agg * phi + k_agg2 * 3 * phi**2
    aggregation_potential = KBT * agg_potential
    mechanical_potential = mismatch_potential + entropy_potential + aggregation_potential
    return mechanical_potential    

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

def calc_mu0(mu0_exp):
    mu0RT = mu0_exp / RT_kcal_mol
    mu0 = mu0RT * RT
    return mu0

def chem_pot_partition(phitot, conc, alength , theta, mu0, area):
    rm = alength / theta
    V_dome = area * thickness
    mem_conc = phitot / (Na * V_dome)
    mu = RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    phi = phitot / (area * phisat)
    return (mu + mu0) * phitot  #This should give a value in J/molecule.


def aggregation_energy(ph_en, phitot, area, k_agg, k_agg2):
    polynomial_eval = ph_en + 0.5 * k_agg * ph_en**2 + k_agg2 * ph_en**3
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

def total_energy(x, phitot, mu0_exp, k_agg):
    k_agg2 = -k_agg / 2.0
    alength = x[0] * 1e-6
    theta = x[1]
    rp = x[2] * 1e-6
    area = calc_dome_area(alength, theta)
    ph_en = phitot / (area * phisat)
    mismatch_E = mismatch( alength, theta, ph_en, phitot, area )
    entropy_E = entropy( ph_en, phitot, area, khat)
    bending_E = dome_bending(alength, theta, ph_en, area)
    tension_E = fm_tension( ph_en, phitot, area )
    agg_en = aggregation_energy(ph_en, phitot, area, k_agg, k_agg2)
    mu0 = calc_mu0(mu0_exp)
    conc = (phi_entire - phitot) / ( V_cyt * Na )
    chem_pot = chem_pot_partition(phitot, conc, alength ,  theta, mu0, area)
    
    saddle_tension = fp_tension(alength, theta, rp)
    saddle_bend = Fplus(alength, theta, rp)
    zero_surface = Fzero(alength, theta, rp)
    cross_section = np.pi * (0.5 * dendDia)**2
    pressure_volume_work = pv_work(alength, theta, rp)
    constraint_rad_gap = constraint_yR(alength, theta, rp)
    constraint_area_energy = area_constraint(alength, theta, rp)
    total_energy = mismatch_E + entropy_E + bending_E + tension_E + agg_en + chem_pot + saddle_tension + saddle_bend + constraint_rad_gap
    return total_energy / KBT

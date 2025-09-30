from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
import numpy as np
import dendShape
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#import seaborn as sns
import pandas as pd
import tag
#sns.set_context("poster")

KBT = 1.380649e-23 * 300
Na = 6.022140e23
RT = 2479 / Na
RT_kcal_mol = 0.593 #kcal/mol
k = 60 * KBT
khat = 40 * KBT
kentropy = KBT
thickness = 5e-9
fm = 0.0e-6
dendDia = 0.5 * (tag.bigDia + tag.smallDia) * 1e-6
Rc = 0.5 * dendDia
k_bind_neq = 0.0
k_bind_eq = 0.0
khat_KBT = khat/KBT
phisat =  1/50e-18
fs = 55 * 1e-7
Cp = 1 / 20e-9
f_cm = 1
Length = tag.Length
xtot = (Length) * 1e-6       
diffL = 20e-9
zero_sub = 1.0
#x0 = [0.5, 0.5, 0.5]  #Note that the init. guess of Rp should be <=0.5 for consistency across all initial guesses.
x0 = np.array([0.5, 0.5, 0.5])

rmbounds = (0.001,5.0)
rpbounds = (0.001,5.0)
tbounds = (0.0001, 1.5)
nbounds = (1.0, 10.0)
part_coeff = 1.0
mu0_exp = tag.mu0_exp
k_agg = tag.k_agg
k_agg2 = tag.k_agg2
k_agg3 = 0.0
V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6
V_voxel = np.pi * (0.5 * dendDia)**2 * diffL
#mu0_exp = -8.6 #kcal/mol Takemura Nat. Sci. Rep., 2017
print("RT: ", RT)
print("KBT: ", KBT)

def calc_sagit(rm, theta, rp):
    r0 = (rm + rp) * np.sin(theta)
    eq_coeff = [1, -2 * Rc, 0.25 * (2 * r0)**2]
    sol = np.roots(eq_coeff)
    sagitta = np.real(sol)[1]
    return sagitta

def calc_dome_area(rm, theta):
    return 2 * np.pi * rm**2 * (1 - np.cos(theta))

def calc_saddle_area(rm, theta, rp):
    r0 = (rp + rm) * np.sin(theta)
    return 2 * np.pi * rp * ( r0 * theta + rp * ( np.cos(theta) - 1 ) )

def calc_mu0():
    mu0RT = mu0_exp / RT_kcal_mol
    mu0 = mu0RT * RT
    return mu0

def chem_pot_partition(phitot, conc, rm , mu0, area):
    V_dome = area * thickness
    mem_conc = phitot / (Na * V_dome)
    #print("Area: ", area)
    #print("MEM CONC: ", mem_conc)
    mu = RT * np.log( (mem_conc) / (conc) )
    phi = phitot / (area * phisat)
    return (mu + mu0) * phitot  #This should give a value in J/molecule.


def aggregation_energy(ph_en, phitot, area):
    polynomial_eval = ph_en + 0.5 * k_agg * ph_en**2 + k_agg2 * ph_en**3
    #polynomial_eval = 0.5 * k_agg * ph_en**2 
    return area * (kentropy * phisat) * polynomial_eval


def entropy( ph_en, phitot, area, khat):
    entropy = area * (kentropy * phisat) * (ph_en * np.log(ph_en) + (1 - ph_en) * np.log(1 - ph_en))
    return entropy

def mismatch( rm, phi, phitot, area ):
    # Note that the  area term cancels out with the phitot/(area)
    #return area * 0.5 * khat * phi * (2/rm - Cp) * (2/rm - Cp)
    return 0.5 * khat * (phitot / phisat) * (1/rm - Cp) * (1/rm - Cp)

def fm_tension( ph_en, phitot, area ):
    return area * fs

def fp_tension(rm, theta, rp):
    saddle_area = calc_saddle_area(rm, theta, rp)
    return saddle_area * fs

def Fplus(rm, theta, rp):
    r0 = (rp + rm) * np.sin(theta)
    r_psi = rm * np.sin(theta)
    saddle_area = calc_saddle_area(rm, theta, rp)
 

    mean_bend_energy = saddle_area * 0.5 * k * ( ( 1/rp + 1/r_psi ) / 2.0 )**2
    gauss_bend_energy = saddle_area * 0.5 * k * 1 / (rp * r_psi)

    mean_energy = fs * np.pi * r0**2

    #total_bend_energy = mean_bend_energy + gauss_bend_energy - mean_energy
    total_bend_energy = mean_bend_energy - mean_energy + gauss_bend_energy

    return total_bend_energy

def constraint_yR(rm, theta, rp):
    r0 = (rm + rp) * np.sin(theta)
    saddle_area = calc_saddle_area(rm, theta, rp)
    dome_area = calc_dome_area(rm, theta)
    total_area = saddle_area + dome_area
    gamma = np.arctan(r0 / Rc)
    yR = Rc / np.cos(gamma)
    cost_var = (yR - Rc)**2 / Rc**2
    return cost_var * total_area * fs

def Fzero(rm, theta, rp):
    r0 = (rm + rp) * np.sin(theta)
    return fs * np.pi * r0**2


def constraint_abs(rm, theta):
    return fs * np.pi * rm**2 * (1 - np.cos(theta))**2

#def constraint_abs(rm, theta):  #Included twice the tension term. So no need to add the tension energy in the total_energy
#    return fs * np.pi * rm**2 * (1 - np.cos(theta)) * (3 - np.cos(theta))

def constraint_ratio(rm, theta):
    return fs * np.pi * rm**2 * (1 - np.cos(theta))**3 / (1 + np.cos(theta))

def dome_bending(rm, ph_en, area):
    return area * (k / (2 * rm * rm))

def calc_dome_volume(rm, theta):
    sphere_rad = rm * np.sin(theta)
    return (2 * np.pi / 3.0) * sphere_rad**3

def calc_saddle_volume(rm, theta, rp):
    r0 = (rm + rp) * np.sin(theta)
    first_term = r0**2 * rp * (1 - np.cos(theta))
    second_term = r0 * rp**2 * (theta - np.sin(theta) * np.cos(theta))
    third_term = (rp**3 / 12.0) * (np.cos(3 * theta) - 9 * np.cos(theta) + 8)
    volume = np.pi * (first_term + second_term + third_term)
    return volume

def pv_work(rm, theta, rp):
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

def constraint_sagitt(rm , theta, rp):
    saddle_area = calc_saddle_area(rm, theta, rp)
    dome_area = calc_dome_area(rm, theta)
    total_area = saddle_area + dome_area
    deltaH = calc_sagit(rm, theta, rp)
    return (deltaH / dendDia) * fs * total_area

def total_energy(x, phitot, conc):
    alength = x[0] * 1e-6
    theta = x[1]
    rm = alength / theta
    rp = x[2] * 1e-6
    area = calc_dome_area(rm, theta)
    ph_en = phitot / (area * phisat)
    mismatch_E = mismatch( rm, ph_en, phitot, area )
    entropy_E = entropy( ph_en, phitot, area, khat)
    bending_E = dome_bending(rm, ph_en, area)
    tension_E = fm_tension( ph_en, phitot, area )
    agg_en = aggregation_energy(ph_en, phitot, area)
    mu0 = calc_mu0()
    chem_pot = chem_pot_partition(phitot, conc, rm ,  mu0, area)
    
    saddle_tension = fp_tension(rm, theta, rp)
    saddle_bend = Fplus(rm, theta, rp)
    zero_surface = Fzero(rm, theta, rp)
    cross_section = np.pi * (0.5 * dendDia)**2
    pressure_volume_work = pv_work(rm, theta, rp)
    #constraint_energy = constraint_sagitt(rm, theta, rp)
    constraint_rad_gap = constraint_yR(rm, theta, rp)
    total_energy = mismatch_E + entropy_E + bending_E + tension_E + agg_en + chem_pot + saddle_tension + saddle_bend + constraint_rad_gap
    return total_energy / KBT

def mech_potential(hmean, phi):
    mismatch_potential = ( 0.5 * khat / phisat ) * (hmean - Cp)**2
    entropy_potential = KBT * ( np.log( phi ) - np.log( 1 - phi ) )
    agg_potential = 1 + k_agg * phi + k_agg2 * 3 * phi**2
    aggregation_potential = KBT * agg_potential
    mechanical_potential = mismatch_potential + entropy_potential + aggregation_potential
    return mechanical_potential

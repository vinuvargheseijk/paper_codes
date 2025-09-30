import numpy as np

thickness = 5e-9
fs = 55 * 1e-7
phisat =  1/50e-18
KBT = 1.380649e-23 * 300
Na = 6.022140e23
RT = 2479 / Na
RT_kcal_mol = 0.593 #kcal/mol
k = 60 * KBT
khat = 40 * KBT
kentropy = KBT
Cp = 1 / 20e-9
dendDia = 1e-6
Rc = 0.5 * dendDia

def calc_dome_area(rm, theta):
    return 2 * np.pi * rm**2 * (1 - np.cos(theta))

def calc_saddle_area(rm, theta, rp):
    r0 = (rp + rm) * np.sin(theta)
    return 2 * np.pi * rp * ( r0 * theta + rp * ( np.cos(theta) - 1 ) )

def calc_mu0(mu0_exp):
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


def aggregation_energy(ph_en, phitot, area, k_agg):
    k_agg2 = -k_agg/2.0
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

def total_energy(x, phitot, conc, k_agg, mu0_exp):
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
    agg_en = aggregation_energy(ph_en, phitot, area, k_agg)
    mu0 = calc_mu0(mu0_exp)
    chem_pot = chem_pot_partition(phitot, conc, rm ,  mu0, area)
    
    saddle_tension = fp_tension(rm, theta, rp)
    saddle_bend = Fplus(rm, theta, rp)
    zero_surface = Fzero(rm, theta, rp)
    cross_section = np.pi * (0.5 * dendDia)**2
    constraint_rad_gap = constraint_yR(rm, theta, rp)
    total_energy = mismatch_E + entropy_E + bending_E + tension_E + agg_en + chem_pot + saddle_tension + saddle_bend + constraint_rad_gap
    return total_energy / KBT


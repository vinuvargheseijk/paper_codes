import numpy as np

part_coeff = 1
Na = 6.022140e23
RT = 2479 / Na
phisat = 1/50e-18 
RT_kcal_mol = 0.593 #kcal/mol
Length = 20
dendDia = 1e-6
thickness = 5e-9
V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6
Na = 6.022140e23
KBT = 1.380649e-23 * 300
khat = 40 * KBT
Cp = 1 / 20e-9
Kb = 40.0

def calc_mu0(mu0_exp):
    mu0RT = mu0_exp / RT_kcal_mol
    mu0 = mu0RT * RT
    return mu0

def mech_potential_calc(hmean, phi, k_agg, k_agg2):
    mismatch_potential = ( 0.5 * khat / phisat ) * (hmean - Cp)**2
    entropy_potential = KBT * ( np.log( phi ) - np.log( 1 - phi ) )
    agg_potential = 1 + k_agg * phi + k_agg2 * 3 * phi**2
    aggregation_potential = KBT * agg_potential
    mechanical_potential = mismatch_potential + entropy_potential + aggregation_potential
    return mechanical_potential    

def chem_pot_partition(phitot, conc, alength , theta, mu0, area):
    rm = alength / theta
    V_dome = area * thickness
    mem_conc = phitot / (Na * V_dome)
    mu = RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    phi = phitot / (area * phisat)
    return (mu + mu0) * phitot  #This should give a value in J/molecule.

def Kf_calc(phi_entire, phitot, rm, phi, mu0_exp, theta, k_agg):
    mu0 = calc_mu0(mu0_exp)
    hmean = 1/rm
    conc = (phi_entire - phitot) / (V_cyt * Na)
    alength = rm * theta
    area = 2 * np.pi * rm**2 * (1 - np.cos(theta))
    chem_potential = chem_pot_partition(phitot, conc, alength, theta, mu0, area) / phitot  #Potential is per molecule
    mech_potential = mech_potential_calc(hmean, phi, k_agg, -k_agg / 2.0)
    Kf = Kb * part_coeff * np.exp( -(chem_potential + mech_potential) / RT )
    print("Kf: ", Kf)
    print("Chem_potential, Mech_potential: ", chem_potential, mech_potential)
    return Kf, chem_potential, mech_potential

from matplotlib.widgets import Slider, Button
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as mplcm
import matplotlib.colors as colors
from mpl_toolkits import mplot3d
import sys
import rdesigneur as rd
import moose
import argparse
from scipy import optimize
import cProfile
from matplotlib.widgets import TextBox
from colorama import Fore
import os
import os.path
import xml.etree.ElementTree as ET

# IMPORTANT CHECKS : UNIT COVERSION IN EXPONENTIAL (KBT/RT)
#                  : MEMBRANE IS ONLY DISTRIBUTED ONCE. HENCE IN SUBSEQUENT STEPS, WITH DOME CHANGES THERE STILL WILL BE SOME RESIDUAL MOLECULES LEFT IN NON-DOME REGION - INCRE                     ASE THE TURNOVER"
#                  : length_under_dome in chem_pot_partition
plot_interval = 1
KBT =  1.380649e-23 * 300
Na = 6.022140e23
RT = 2479 / Na  #J/molecule
k = 63 * KBT
khat = 35 * KBT    # 35+- 7 k_BT       35 * 1.380649e-23 * 300    Joules
kentropy = 1 * KBT # The full term is kentropy * kBT*phitot * R^2
thickness = 10e-9
dist_mem = 0.5e-6
dist_near_mem = 0.45e-6
dist_recruit = dist_mem - dist_near_mem
fm = 0.0e-6
dendDia = 1e-6
k_bind_neq = 0.0
k_bind_eq = 0.0
k_agg = -20
khat_KBT = khat/KBT
phisat = 1/50e-18
fs = 55 * 1e-7
diffL = 10e-9
diffDt = 0.0005
Cp = 1/20e-9
f_cm = 1
p_start = 5
p_end = 15
num_lines = 50
Length = 15
xtot = (Length) * 1e-6       
constr_multiplier = 1.0
plot_color = ['b', 'g', 'r', 'c', 'm', 'y', 'k'] 
disable_gui = True
x0 = [1.0, 0.5]
rmbounds = (0.001,10.0)
tbounds = (0.001,3.14)
nbounds = (1.0, 10.0)
#mu0 = -10 * KBT
p = 1
Lat_diff = 1.0e-12
V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6

#OLD conversion
RT_kJpM = 2.479  #RT in kJ / mol
RT_Jpmolecule = RT_kJpM * 1e3 / Na  #RT in J/molecule
Jpmolecule = 1 / RT_Jpmolecule  #inverse of previous. Since RT is in the denominator, this goes into numerator

#NEW conversion
mu0 = -8.6 #kcal/mol Takemura Nat. Sci. Rep., 2017
mu0RT = mu0 * (1 / 0.593) #1 kcal/mol is 1/0.593 RT. Now the unit is in RT
mu0 = mu0RT * RT #in J/molecule




Lat_diff_time = (0.5 * dendDia)**2 / (2 * Lat_diff)
Kb = 1/Lat_diff_time
print("Kb: ", Kb)

#Length = Length + space_val
#phisat = 50e-18 # Saturation level: One molecule in 50 square nanometers

def writeXML_td(time_d,fileName_td):
    if os.path.isfile(fileName_td):
       tree = ET.parse(fileName_td)
       root = tree.getroot()
       #distr = ET.SubElement(root, 'calcium'+str(sys.argv[3])+"oth"+str(sys.argv[4]))
       for i in range(len(time_d)):
          avec = ET.SubElement(root, 'time_d'+str(i)+str(fileName_td))
          avec.text = ''.join(str(j)+' ' for j in time_d[i])+'\n'
       tree.write(fileName_td)
    else:
       root = ET.Element('Data')
       for i in range(len(time_d)):
          avec = ET.SubElement(root, 'time_d'+str(i)+str(fileName_td))
          avec.text = ''.join(str(j)+' ' for j in time_d[i])+'\n'
       #times = ET.SubElement(root, 'times')
       #times.text = ''.join(str(j)+' ' for j in timeArr)+'\n'
       #xmaxFWHH = ET.SubElement(root, 'xmaxFWHH'+trialNum)
       #xmaxFWHH.text = ''.join(str(k)+' ' for k in maxFWHH)+'\n'
       tree = ET.ElementTree(root)
       tree.write(fileName_td)

def mid_diff_profile(x,t,diff, strength):
    print("Strength", strength)
    y = [0] * len(x)
    print("Length of x", len(x))
    print("Length of y", len(y))
    for i in range(len(x)):
        y[i] = (strength/(t)**0.5) * np.exp(-x[i]**2/(4 * diff * t))
    return y

def calc_saddle_area(rm, rp, theta):
    r0 = (rp + rm) * np.sin(theta)
    return 2 * np.pi * rp * (r0*theta + rp * ( np.cos( theta ) - 1 ) )

def entropy( area, phitot, phisat, khat, fs, k ):
    ph_en = phitot/ (area * phisat)
    entropy = (kentropy * phisat) * (ph_en*np.log(ph_en) - ph_en)
    return entropy

def aggregation_energy(rm, theta, phitot, phisat, k_agg_v, area):
    ph_en = phitot/ (phisat * area)
    return area * (kentropy * phisat) * ( 0.5 * k_agg_v * ph_en**2 + 5 * ph_en**3 + 6.12 * ph_en**4 + 7.06 * ph_en**5)



def mismatch( rm, phitot, phisat, khat, fs, k ):
    # Note that the  area term cancels out with the phitot/(area)
    return 0.5 * khat * (phitot / phisat) * (1/rm - Cp) * (1/rm - Cp)

def bind_energy(phitot, phisat, area):
    ph_en = phitot/ (phisat * area)
    return area * KBT * k_bind_neq * phisat * ph_en

def fp_tension( rp, rm, theta, khat, fs, k ):
    r0 = (rp + rm) * np.sin( theta )
    #return 2 * np.pi * rp * fs * (r0*theta + rp * ( np.cos( theta ) - 1 ) )
    return calc_saddle_area(rm, rp, theta) * fs

def fm_tension( rm, theta, phitot, phisat, khat, fs, k ):
    area = 2 * np.pi * rm * rm * ( 1 - np.cos( theta ) )
    return area * fs

def Fplus( rp, rm, theta, khat, fs, k ):
    r0 = (rp + rm) * np.sin( theta )
    r_psi = r0 - rp * np.sin(theta) 
    mean_curv = (1/r_psi + 1/rp)/2.0
    area_s = calc_saddle_area(rm, rp, theta)
    mean_curv_energy = area_s * 0.5 * k * mean_curv**2
    ret = 2 * np.pi * (k * theta ) + mean_curv_energy
    return ret

def Fminus( rm, theta, phitot, phisat, khat, fs, k ):
    area = 2 * np.pi * rm * rm * ( 1 - np.cos( theta ) )
    Fm = mismatch( rm, phitot, phisat, khat, fs, k) + area * ( k/(2*rm*rm) + entropy( area, phitot, phisat, khat, fs, k ) )
    return Fm

def Fzero( rp, rm, theta, n, khat, fs, k, order_spine, pXbend):
    excess_area, total_spine_area, original_area = membrane_difference_area(rp, rm, theta)
    return excess_area * fs

def constraint(rp , rm, theta, dendDia, fs, phitot, area):
    ph_en = phitot/ (phisat * area)
    return fs *  ( (rp + rm) * (1 - np.cos(theta)) ) * ( (rp + rm) * (1 - np.cos(theta)) + (dendDia/2.0) ) * ( (rp + rm) * (1 - np.cos(theta)) + (dendDia/2.0) )/dendDia

def area_constraint(rp , rm, theta, dendDia, fs, phitot, area):
    return fs * area * rm/dendDia

def membrane_difference_area(rp, rm, theta):
    r0 = (rp + rm) * np.sin( theta )
    dome_area = 2 * np.pi * rm * rm * (1 - np.cos(theta))
    saddle_area = calc_saddle_area(rm, rp, theta)
    total_spine_area = dome_area + saddle_area
    original_area = np.pi * r0 * r0
    return (total_spine_area - original_area), total_spine_area, original_area




def total_energy(x,phitot, phisat, khat, fs, k, n, dummy, order_spine, pXbend, conc):
    n = 1
    mu = chem_pot_partition(phitot, p, conc, x[0], x[1])
    rm = x[0] * 1e-6
    theta = x[1]
    rp = rm * np.sin(theta)
    phitot = phitot/n
    area = 2 * np.pi * rm * rm * ( 1 - np.cos( theta ) )
    phi = phitot / ( area * phisat )
    Fp = Fplus( rp, rm, theta, khat, fs, k )
    Fm = Fminus( rm, theta, phitot, phisat, khat, fs, k )
    F0 = Fzero( rp, rm, theta, n, khat, fs, k, order_spine, pXbend )
    try:
       k_agg_v = x[2]
    except:
       k_agg_v = k_agg 
    agg_en = aggregation_energy(rm, theta, phitot, phisat, k_agg_v, area)


    basal_tension_energy = ( fs + fm ) * 2 * np.pi * dendDia * xtot/int(num_spines_v.text_disp._text)

    total = n * ( Fm + Fp ) + membrane_tension(rm,theta,phisat,phitot,n,khat,fs, k, order_spine, pXbend) + n * agg_en + constr_multiplier * n * area_constraint(rp , rm, theta, dendDia, fs, phitot, area) + mu * phitot

    #Check if constraint is included
    return total * 1e16

def chem_pot_partition(phitot, p, conc, rm , theta):
    area_dome = 2 * np.pi * (rm * 1e-6)**2 * (1 - np.cos(theta))
    V_dome = area_dome * thickness
    mem_conc = phitot / (Na * V_dome)
    mu = RT * np.log( (mem_conc + 1e-24) / (conc + 1e-24) )
    return (mu + mu0)

def chem_pot_calc(memC_zone, cytC_zone):
    mu = mu0 + RT * np.log( (memC_zone + 1e-24) / (cytC_zone + 1e-24) )
    return mu


def full_potential(phitot, rm, theta, conc, memConc):
    area = 2 * np.pi * (rm * 1e6)**2 * (1 - np.cos(theta))
    phi = phitot / (area * phisat)
    mu_nernst = chem_pot_partition(phitot, 1, conc, rm , theta)
    #mu_nernst = chem_pot_calc(memConc, conc)
    mu = (0.5 * khat / phisat) * (1 / rm - Cp)**2 + KBT * np.log(phi) + KBT * (phi + 15 * phi**2 + 24.48 * phi**3 + 35.3 * phi**4) + mu_nernst
    return mu

def membrane_tension(x,y,phitot,phisat,n,khat,fs, k, order_spine, pXbend):
    rm = x
    theta = y
    rp = rm * np.sin(theta)
    F0 = Fzero( rp, rm, theta, n, khat, fs, k, order_spine, pXbend )
    fm_ten = n * fm_tension( rm, theta, phitot, phisat, khat, fs, k )
    fp_ten = n * fp_tension( rp, rm, theta, khat, fs, k )
    mem_ten = F0 + fm_ten + fp_ten
    return mem_ten

def bending(x,y,phitot,n,khat,fs):
    rm = np.asarray(x)
    theta = y
    rp = rm * np.sin(theta)
    r0 = (rp + rm) * np.sin( theta )
    r_psi = r0 - rp * np.sin(theta)
    mean_curv = (1/r_psi + 1/rp)/2.0
    area_s = calc_saddle_area(rm, rp, theta)
    mean_curv_energy = area_s * 0.5 * k * mean_curv**2

    area_m = 2 * np.pi * rm * rm * ( 1 - np.cos( theta ) )
    fp_bend = 2 * np.pi * k * theta + mean_curv_energy
    fm_bend = area_m *  k/(2*rm*rm)
    mem_bend = fp_bend + fm_bend
    return mem_bend


def chem_model(khat, k_agg, k, phisat, steady):
    moose.Neutral( '/library' )
    moose.Neutral( '/library/cm' )
    moose.CubeMesh( '/library/cm/dend' )
    cytosol = moose.Pool( '/library/cm/dend/cytosol' )
    cytConc = moose.Pool( '/library/cm/dend/cytConc' )
    memConc = moose.Pool( '/library/cm/dend/memConc' )
    memv = moose.Pool( '/library/cm/dend/memv' )
    membrane = moose.Pool( '/library/cm/dend/membrane' )
    Hmean = moose.Pool( '/library/cm/dend/Hmean' )
    Sphi = moose.Pool( '/library/cm/dend/Sphi' )
    marker = moose.Pool( '/library/cm/dend/marker' )
    theta = moose.Pool( '/library/cm/dend/theta' )
    pot_eval = moose.Pool( '/library/cm/dend/pot_eval' )
    pot_nernst = moose.Pool( '/library/cm/dend/pot_nernst' )


    reacF = moose.Reac( '/library/cm/dend/reacF' )
    reacB = moose.Reac( '/library/cm/dend/reacB' )

    cytosol.diffConst = 1e-12
    membrane.diffConst = 0.0   #Check
    cytosol.concInit = 0.0

    moose.connect( reacF, 'sub', cytosol, 'reac')
    moose.connect( reacF, 'prd', membrane, 'reac')
    moose.connect( reacB, 'sub', membrane, 'reac')
    moose.connect( reacB, 'prd', cytosol, 'reac')

    reacF.Kf = 0
    reacF.Kb = 0
    reacB.Kf = 0
    reacB.Kb = 0
    func_kf = moose.Function('/library/cm/dend/func_kf')
    func_pot_eval = moose.Function('/library/cm/dend/func_pot_eval')
    func_pot_nernst_eval = moose.Function('/library/cm/dend/func_pot_nernst_eval')
    func_pot_eval.x.num = 4
    func_pot_nernst_eval.x.num = 3
    func_kf.x.num = 2

    voxel_volume = np.pi * (0.5 * dendDia)**2 * diffL

    func_kf.expr = str(Kb) + " * x0 * exp( -1 * (x1 / " + str(KBT) + ") )"  


    func_pot_eval.expr = "x0 * (" + str( 0.5 * khat / phisat ) + " * (x1 - " + str(Cp) + ")^2 + " + str(KBT) + " * log(x2 + 1e-2) + " + str(KBT) + " * (x2 + 15 * x2^2 + 24.48 * x2^3 + 35.3 * x2^4) + x3 )" 

    func_pot_nernst_eval.expr = "x0 * ( " + str(mu0) + " + " + str(RT) + " * log( (x1 + 1e-24) / (x2 + 1e-24) ) )"


    moose.connect( marker, 'nOut', func_kf.x[0], 'input' ) 
    moose.connect( pot_eval, 'nOut', func_kf.x[1], 'input' )
    moose.connect( func_kf, 'valueOut', reacF, 'setNumKf' )
    
    moose.connect( marker, 'nOut', func_pot_eval.x[0], 'input' ) 
    moose.connect( Hmean, 'nOut', func_pot_eval.x[1], 'input' )
    moose.connect( Sphi, 'nOut', func_pot_eval.x[2], 'input' )
    moose.connect( pot_nernst, 'nOut', func_pot_eval.x[3], 'input' )
    moose.connect( func_pot_eval, 'valueOut', pot_eval, 'setN' )

    moose.connect( marker, 'nOut', func_pot_nernst_eval.x[0], 'input' ) 
    moose.connect( memConc, 'nOut', func_pot_nernst_eval.x[1], 'input' ) 
    moose.connect( cytConc, 'nOut', func_pot_nernst_eval.x[2], 'input' ) 
    moose.connect( func_pot_nernst_eval, 'valueOut', pot_nernst, 'setN' )

   
   
    reacB.Kf = Kb

    reacB.Kb = 0.0

    reacF.Kb = 0.0

    
 


def rdes_build(steady):
    #cm = chem_model(khat, k_agg, k, phisat,Hmean_val, Hmean_val2, f_m_val, f_m2_val, phitot, theta_val, theta_val2, steady, phi_1, phi_2)
    cm = chem_model(khat, k_agg, k, phisat, steady)
    rdes = rd.rdesigneur(
    turnOffElec = True,
    diffusionLength = diffL,
    diffDt = diffDt,
    chemDt = diffDt,
    #Check Length
    #cellProto = [['ballAndStick', 'soma', 1e-6, 0e-6, dendDia, Length * 1e-6, 1]],
    cellProto = [['somaProto', 'soma', dendDia, (Length) * 1e-6]],  #Check
    chemProto = [['cm', 'cm']],
    chemDistrib = [['cm', 'soma', 'install', '1' ]],    
    #chemDistrib = [['cm', 'dend#', 'install', '1' ]],    
    )
    print("Building model")
    model = rdes.buildModel('/model')
    print(moose.le('/model/elec'))
    moose.setClock(30, diffDt)   #PyRun clock tick 30
    moose.setClock(10, diffDt)
    for i in range( 11, 18 ):
       moose.setClock( i, diffDt )
    moose.setClock( 18, diffDt )

def find_shape(x0, rmbounds, tbounds, nbounds, phitot, num_spines, n, order_spine, pXbend, conc):
    ret = optimize.minimize( total_energy, x0=x0, bounds = (rmbounds,tbounds),args=(phitot,phisat, khat, fs, k, n, 0, order_spine, pXbend, conc),method = 'L-BFGS-B', options = {'ftol':1e-12} )
    theta = ret.x[1]
    rm = ret.x[0]
    Hmean_val = 1/(ret.x[0] * 1e-6)
    return rm, theta, Hmean_val, n

def corrupt_opt_check(rm, theta):
    break_flag = 0
    if np.abs(rm - rmbounds[0]) < 1e-6 or np.abs(rm - rmbounds[-1]) < 1e-6 or np.abs(theta - tbounds[0]) < 1e-6 or np.abs(theta - tbounds[-1]) < 1e-6:
        break_flag = 1
    #if np.abs(rm - x0[0]) < 1e-6 or np.abs(theta - x0[1]) < 1e-6:
    #    break_flag = 1
    return break_flag    

def distributing_shape(x0, rmbounds, tbounds, nbounds, phitots, num_spines, spacing, x2, y2, n, x_positions_initial, run_time, concs):
    dome_starts = []
    dome_ends = []
    shape_rm = []
    shape_theta = []
    shape_rp = []
    shape_phi = []
    prev_xbend = []
    prev_Lflat = (Length - ( (int(num_spines_v.text_disp._text) - 1) * float(space_v.text_disp._text) ) )/2.0
    prev_xbend.append(0)
    x_positions = []
    for shape in range(int(num_spines_v.text_disp._text)):
        print("CONC: ", concs[shape])
        rm, theta, Hmean_val, n = find_shape(x0, rmbounds, tbounds, nbounds, phitots[shape], int(num_spines_v.text_disp._text),n, shape, prev_xbend[shape], concs[shape])
        break_flag = corrupt_opt_check(rm, theta)
        if break_flag == 1:
            if disable_gui == True:
              print("Optimization didn't work for phitot: ", phitots[shape])
              print("RM, THETA: ", rm, theta)
            if disable_gui == False:
              print(Fore.RED + "Optimization didn't work for phitot: ", phitots[shape])
              print("RM, THETA: ", rm, theta)
            break
        rp = rm * np.sin(theta)
        new_dome_area = 2 * np.pi * (rm * 1e-6)**2 * (1 - np.cos(theta))
        f_m_val = new_dome_area
        phi = phitots[shape]/(f_m_val * phisat)
        shape_rm.append(rm)
        shape_rp.append(rp)
        shape_theta.append(theta)
        shape_phi.append(phi)
        xbend = 2 * np.sin( theta ) * (rp + rm )
        prev_xbend.append(xbend * 1e-6)
        x, y, dome_start, dome_end, Lflat, x_pos = dendShape( theta, rp, rm, shape, prev_Lflat, x_positions_initial, run_time )
        print("Position of spine" +str(shape) +":"+str(x_pos))
        x_positions.append(x_pos)
        prev_Lflat = Lflat + float(space_v.text_disp._text)
        dome_starts.append(dome_start)
        dome_ends.append(dome_end)
        x2 = np.append(x2, x)
        y2 = np.append(y2, y)
    return x2, y2, f_m_val, dome_starts, dome_ends, new_dome_area, phi, shape_rm, shape_rp, shape_theta, shape_phi, prev_xbend, x_positions 

def marking(y):
    print("MARKING: ", y)
    marking.action = 'False'


def finding_dome(num_voxels, phitots, dome_starts, dome_ends):
       marker = [0] * num_voxels
       for i in range(num_voxels):
           marker[i]=0
           for j in range(int(num_spines_v.text_disp._text)):
              if (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) >= dome_starts[j] and (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) <= dome_ends[j]:
                     marker[i] = 10e-1
       moose.vec( '/model/chem/dend/marker' ).n = marker          
       return marker

def allot_geom_par(num_voxels, dome_starts, dome_ends, shape_rm, shape_rp, shape_theta, shape_phi, marker):
       thetaVec = moose.vec( '/model/chem/dend/theta' ).n
       HmeanVec = moose.vec( '/model/chem/dend/Hmean' ).n
       SphiVec = moose.vec( '/model/chem/dend/Sphi' ).n
       for i in range(num_voxels):
           for j in range(int(num_spines_v.text_disp._text)):
               if (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) >= dome_starts[j] and (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) <= dome_ends[j]:
                    thetaVec[i] = marker[i] * shape_theta[j]
                    HmeanVec[i] = marker[i] * 1/(shape_rm[j] * 1e-6)
                    SphiVec[i] = marker[i] * shape_phi[j]
       moose.vec( '/model/chem/dend/theta' ).n = thetaVec
       moose.vec( '/model/chem/dend/Hmean' ).n = HmeanVec
       moose.vec( '/model/chem/dend/Sphi' ).n = SphiVec

def allot_conc(num_voxels, marker, C_zones, dome_starts, dome_ends, name):
    Vec = moose.vec( '/model/chem/dend/' + str(name) ).n
    for i in range(num_voxels):
        for j in range(int(num_spines_v.text_disp._text)):
            if (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) >= dome_starts[j] and (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) <= dome_ends[j]:
                Vec[i] = marker[i] * C_zones[j]
    moose.vec( '/model/chem/dend/' +str(name) ).n = Vec             

def  allot_volume(num_voxels, marker, dome_starts, dome_ends, name, spine_length_distribution, dome_volume_distribution):
    Vec = moose.vec( '/model/chem/dend/' + str(name) ).n
    for i in range(num_voxels):
        for j in range(int(num_spines_v.text_disp._text)):
            if (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) >= dome_starts[j] and (moose.element('/model/chem/dend/mesh['+str(i)+']').Coordinates[0]) <= dome_ends[j]:
                
                Vec[i] = marker[i] * (dome_volume_distribution[j] / spine_length_distribution[j])
    moose.vec( '/model/chem/dend/' +str(name) ).n = Vec             




def spine_content(membrane, concs, dome_starts, dome_ends, mem_coords):
       length_zones = [0] * int(num_spines_v.text_disp._text)
       dome_zones = [0] * int(num_spines_v.text_disp._text)
       memC_zones = [0] * int(num_spines_v.text_disp._text)
       spine_count = 0
       num_voxels = len(membrane)
       
       mem_array = np.array(membrane)
       #mem_coords = np.arange(0, Length * 1e-6, diffL)
       
       mem_non_zero_list = list(np.where(mem_array > 1e-24)[0])
       mem_count = 0
       mem_non_zero_list.append(0)  #padding with zero so that index will not be run over
       print("LENGTH OF NON ZERO MEM: ", len(mem_non_zero_list))
       for i in range(len(mem_non_zero_list) - 1):
           if mem_non_zero_list[i] + 1 == mem_non_zero_list[i + 1]:
              mem_count = mem_count + membrane[mem_non_zero_list[i]]
           else:
              length_zones[spine_count] = mem_count + membrane[mem_non_zero_list[i]]
              spine_count = spine_count + 1
              mem_count = 0

       for nv in range(num_voxels):
           for ns in range(int(num_spines_v.text_disp._text)):
               if mem_coords[nv] >= dome_starts[ns] and mem_coords[nv] <= dome_ends[ns]:
                   dome_zones[ns] = dome_zones[ns] + membrane[nv]

       print("SPINE COUNT, LENGTH ZONES: ", spine_count, length_zones)            
       print("SUM CYTOSOL: ", sum(moose.vec('/model/chem/dend/cytosol').n))
       print("SUM MEMBRANE: ", sum(moose.vec('/model/chem/dend/membrane').n))
       print("SUM LENGTH ZONES: ", sum(length_zones))
       dome_volume_distro = []
       for pt in range(len(dome_zones)):            
           rm, theta, Hmean_val, n = find_shape(x0, rmbounds, tbounds, nbounds, dome_zones[pt], int(num_spines_v.text_disp._text), 1, 0, 0, concs[pt])       
           dome_volume = 2 * np.pi * (rm * 1e-6)**2 * (1 - np.cos(theta)) * thickness
           memC_zones[pt] = dome_zones[pt] / (Na * dome_volume) 
           dome_volume_distro.append(dome_volume)

       return length_zones, spine_count, memC_zones, dome_volume_distro, dome_zones

def cytosol_content(cytosol, membrane, dome_starts, dome_ends, mem_coords):
    cytN_zones = [0] * int(num_spines_v.text_disp._text)
    cytC_zones = [0] * int(num_spines_v.text_disp._text)
    spine_count = 0
    num_voxels = len(membrane)

    spine_length_distro = []

    for ns in range(int(num_spines_v.text_disp._text)):
        spine_length_distro.append(int( (dome_ends[ns] - dome_starts[ns])/diffL ) )

    """
    mem_array = np.array(membrane)
    mem_non_zero_list = list(np.where(mem_array > 1e-24)[0])
    mem_count = 0
    mem_non_zero_list.append(0)  #padding with zero so that index will not be run over
    
    spine_length = 0
    cytNi = 0
    count = 0
    for i in range(len(mem_non_zero_list) - 1):
           if mem_non_zero_list[i] + 1 == mem_non_zero_list[i + 1]:
              cytNi = cytNi + cytosol[mem_non_zero_list[i]]
              count = count + 1
           else:
              cytN_zones[spine_count] = cytNi + cytosol[mem_non_zero_list[i]]
              cytC_zones[spine_count] = cytNi / (Na * np.pi * (0.5 * dendDia)**2 * spine_length_distro[spine_count] * diffL)
              spine_count = spine_count + 1
              cytNi = 0
              print("COUNTED: ", count)
              print(cytN_zones)
              count = 0
    """

    for ns in range(int(num_spines_v.text_disp._text)):
        cytN_zones[ns] = 0
        cytC_zones[ns] = 0


    for nv in range(num_voxels):
           for ns in range(int(num_spines_v.text_disp._text)):
               if mem_coords[nv] >= dome_starts[ns] and mem_coords[nv] <= dome_ends[ns]:
                   cytN_zones[ns] = cytN_zones[ns] + cytosol[nv]

    for ns in range(int(num_spines_v.text_disp._text)):
           cytC_zones[ns] = cytN_zones[ns] / (Na * np.pi * (0.5 * dendDia)**2 * spine_length_distro[ns] * diffL)
          
    
    return cytN_zones, cytC_zones, spine_length_distro            

   

def pool_homogenize(marker, num_voxels, zones, dome_starts, dome_ends, pool_name, pool_content, nonzero_e1):
       spine_lengths = []
       counting_voxels = 0
       for non in range(len(nonzero_e1[0]) - 1):
                if (nonzero_e1[0][non+1] - nonzero_e1[0][non]) == 1:
                     counting_voxels = counting_voxels + 1
                     print("Counting", counting_voxels)
                else:
                     spine_lengths.append(counting_voxels + 1) # +1 because I didnt count the last voxels of the spine
                     counting_voxels = 0
       spine_lengths.append(counting_voxels + 1) # +1 because I didnt count the last voxel   

       print("SPINE LENGTHS: ", spine_lengths) 
       non_zero = list(nonzero_e1[0])
       non_zero.append(0)
       start = 0
       for length in range(len(spine_lengths)):
           for non in range(start, len(non_zero)):
               print(non)
               if non_zero[non + 1] == non_zero[non] + 1:
                   print(non_zero[non])
                   pool_content[non_zero[non]] = zones[length] / spine_lengths[length]
               else:
                   pool_content[non_zero[non]] = zones[length] / spine_lengths[length]
                   print("BREAK POINT: ", non_zero[non])
                   start = non + 1
                   break
            
       moose.vec('/model/chem/dend/' + str(pool_name)).n = pool_content
       print("SUM OF LENGTH ZONES: ", sum(zones))
       print("SUM MEMBRANE: ", sum(moose.vec('/model/chem/dend/' + str(pool_name)).n))
       print("SPINE LENGTH: ", spine_lengths)
       print("NONZERO: ", len(nonzero_e1[0]))
       print("POOL: ", len(np.where(pool_content > 0)[0]))

       """ 
       
       lengths = []
       counting_voxels = 0
       for non in range(len(nonzero_e1[0]) - 1):
                if (nonzero_e1[0][non+1] - nonzero_e1[0][non]) == 1:
                     counting_voxels = counting_voxels + 1
                     print("Counting", counting_voxels)
                else:
                     lengths.append(counting_voxels + 1) # +1 because I didnt count the last voxels of the spine
                     counting_voxels = 0
       lengths.append(counting_voxels + 1) # +1 because I didnt count the last voxel   
       print("LENGTHS FROM NONZERO: ", lengths)

       
       #Distributing membrane molecules
       start_marker = 0
       for l in range(len(lengths)):
         for non in range(start_marker, start_marker + lengths[l]):
             pool_content[nonzero_e1[0][non]] = dome_zones[l]/lengths[l]
         start_marker = start_marker + lengths[l]
         print("Resetting start marker", start_marker)
         print("Lengths of each spines", lengths)           
         print("POOL NAME: ", pool_name)
       """


def update(val):
    steady = "unsteady_protr"
    rdes_build(steady)
    moose.reinit()

    
   

    scales = list(dome_v.text_disp._text.split(" "))

    for s in range(len(scales)):
        print("SCALES: ", scales)
        scales[s] = float(scales[s])

    num_voxels = len(moose.vec("/model/chem/dend/cytosol").n)
    mem_coords = []
    for nv in range(num_voxels):
        mem_coords.append(moose.element('/model/chem/dend/mesh['+str(nv)+']').Coordinates[0])



    x_diff = np.arange(0,num_voxels,1)
    cyt_diff = moose.element("/model/chem/dend/cytosol").diffConst
    coord_diff = []
    middle_pos = int((0.5 * Length * 1e-6)/diffL)
    coord_diff = (x_diff - middle_pos) * diffL
    stim_center = False
    if stim_center == True:
       y_diff = mid_diff_profile(coord_diff,0.1,cyt_diff, 100)
       phitot = ( int(strength_v.text_disp._text) - sum(y_diff) )/int(num_spines_v.text_disp._text) 
       moose.vec("/model/chem/dend/cytosol").n = y_diff
    else:
       cyt_conc = 0.2e-3
       Vcyt = np.pi * (dendDia/2.0)**2 * Length * 1e-6
       cyt_num = Vcyt * Na * cyt_conc
       phitot = int(strength_v.text_disp._text) - cyt_num 
       initial_phitot = phitot
       #moose.vec("/model/chem/dend/cytosol").n = cyt_num/len(moose.vec("/model/chem/dend/cytosol").n)



        

    #Updating geometry by minimizing for phitot
    dome_starts = []
    dome_ends = []
    dome_start = 0
    done_end = 0
    x = []
    y = []
    x2 = np.asarray(x)
    y2 = np.asarray(y)
    spacing = float(space_v.text_disp._text)
    
    dome_zones = [0] * int(num_spines_v.text_disp._text) 
    for s in range(int(num_spines_v.text_disp._text)): 
       dome_zones[s] = scales[s] * phitot
    print("Scales: ", scales)
    print("PHITOT: ", dome_zones)
    if abs(sum(dome_zones) - phitot) > 0.01:
       print("sum phitot exceeds allotted. Correct the scalings of dome")
       if int(num_spines_v.text_disp._text) != 1:
          if scales[0] !=1: 
             exit()
    if int(num_spines_v.text_disp._text) == 1:
       balance = (1 - scales[0]) * phitot
       print("BALANCE CALC: ", balance, cyt_num, scales[0] * phitot, balance + cyt_num + scales[0] * phitot)
    else:
       balance = 0
    moose.vec("/model/chem/dend/cytosol").n = ( cyt_num + balance ) / len(moose.vec("/model/chem/dend/cytosol").n)
    



    #Finding the initial geometrical parameters. Note that initial phitots are same for all protrusions
    num_voxels = len(moose.vec("/model/chem/dend/cytosol").n)
    concs = [cyt_conc] * int(num_spines_v.text_disp._text)
    x_positions = []
    run_time = 0
    x2, y2, f_m_val, dome_starts, dome_ends, new_dome_area, phi, shape_rm, shape_rp, shape_theta, shape_phi, prev_xbend, x_positions_initial  = distributing_shape(x0, rmbounds, tbounds, nbounds, dome_zones, int(num_spines_v.text_disp._text), spacing, x2, y2, 1, x_positions, run_time, concs)

    dt = 0.025
    prev_X = x2
    prev_Y = y2



    #Check 
    #moose.vec("/model/chem/dend/cytosol").n = ( int(strength_v.text_disp._text) - int(num_spines_v.text_disp._text) * phitot ) / len(moose.vec("/model/chem/dend/cytosol").n)



    runtime_txt = float(runtime_v.text_disp._text)
    k_agg1= k_agg
    k_agg2 = k_agg
                       
    time_list = []
    for i in range(int(num_spines_v.text_disp._text)):
      print(i)  
      globals()["area_" + str(i+1) + "_list"] = []
      globals()["phi_" + str(i+1) + "_list"] = []
      globals()["phitot" + str(i+1) + "_list"] = []
      globals()["rd" + str(i+1) + "_list"] = []
      globals()["mu" + str(i+1) + "_list"] = []
      globals()["Slength" + str(i+1) + "_list"] = []
      globals()["internKf" + str(i+1) + "_list"] = []

    cyt_save = []
    x_save = []
    y_save = []

    if disable_gui == False:
      ax2.set_xlabel("Time (s)")
      ax2.set_ylabel("(Total - Dome) molecules \n in the protrusion")

    for r in range(int((runtime_txt/dt))):
       print("Intial positions are: " + str(x_positions_initial)) 
       currentSimTime = moose.element('/clock').currentTime
       run_time =  currentSimTime
       print("Current Time: ", run_time)
  


       #Finding the voxels under the domes. I assign a variable marker which holds the voxels under the domes. This is assigned to the pool marker in the following function. The pool marker multplies the functions in the chem model. This make sure that that recruitment occur only under the dome region.

       marker = finding_dome(num_voxels, dome_zones, dome_starts, dome_ends)
       nonzero_e1 = np.where(np.array(marker) != 0)                  
       membrane = moose.vec('/model/chem/dend/membrane').n

       #Only at time = 0. I distribute the initial phitot equally for all domes.
       cytosol = moose.vec('/model/chem/dend/cytosol').n    
       if r == 0:
          pool_homogenize(marker, num_voxels, dome_zones, dome_starts, dome_ends, 'membrane', membrane, nonzero_e1)
       """
       """
       
       #Finding the cytosol concentration under the dome
       cytN_zones, cytC_zones, spine_length_distro = cytosol_content(cytosol, membrane, dome_starts, dome_ends, mem_coords)
       print("CYT ZONES: (N, C) ", cytN_zones, cytC_zones)
       print("Spine length distribution: ", spine_length_distro)
       allot_conc(num_voxels, marker, cytC_zones, dome_starts, dome_ends, 'cytConc')
       
       #Now I need to calculate the updated number of molecules under each dome. This is stored in dome_zones.
       length_zones, spine_count, memC_zones, dome_volume_distro, dome_zones = spine_content(membrane, cytC_zones, dome_starts, dome_ends, mem_coords)
           

       if disable_gui == False:  
          ax5.clear() 
          ax10.clear()
          x_vec = np.linspace(0,Length,num_voxels)  #Check
          ax2.plot(r * dt, length_zones[0] - dome_zones[0], "*", color = 'b')    
          if int(num_spines_v.text_disp._text) > 1:
            ax2.plot(r * dt, length_zones[1] - dome_zones[1], "*", color = 'r')    
          ks = moose.element( '/model/chem/dend/ksolve' )
          rates2 = ks.rateVec['/model/chem/dend/reacF']
          ax5.plot(x_vec, rates2, 'r')
          ax5.set_title("$k_{f}$")
          ax6.set_title("Cyt N")
          ax9.set_title("Mem CONC")
          ax6.plot(r * dt, cytN_zones[0], marker = "o", color = 'b')
          if int(num_spines_v.text_disp._text) > 1:
            ax6.plot(r * dt, cytN_zones[1], marker = "X", markersize = 12, markerfacecolor = 'None', color ='r')
          ax9.plot(r * dt, memC_zones[0], marker = "o", color = 'b')
          if int(num_spines_v.text_disp._text) > 1:
            ax9.plot(r * dt, memC_zones[1], marker = "X", markersize = 12, markerfacecolor = 'None', color ='r')
          pot_eval = moose.vec('/model/chem/dend/pot_eval').n
          print("Max of Hmean: ", max(moose.vec('/model/chem/dend/Hmean').n))
          print("Max of Theta: ", max(moose.vec('/model/chem/dend/theta').n))
          print("Max of phi: ", max(moose.vec('/model/chem/dend/Sphi').n))
          ax10.plot(x_vec, pot_eval)
          ax10.set_title("Potential")
          ax.set_title("Kf / Kb in time", fontsize=18)
          print(Fore.GREEN + "Potentials")
          ret = optimize.minimize( total_energy, x0=x0, bounds = (rmbounds,tbounds),args=(dome_zones[0],phisat, khat, fs, k, 1, 0, 0, 0, cytC_zones[0]),method = 'L-BFGS-B', options = {'ftol':1e-12} )
          theta = ret.x[1]
          rm = ret.x[0] * 1e-6
          print( full_potential(dome_zones[0], shape_rm[0], shape_theta[0], cytC_zones[0], memC_zones[0]) )
          print( "nernst: ", chem_pot_calc(memC_zones[0], cytC_zones[0]))
          if int(num_spines_v.text_disp._text) > 1:
            ret = optimize.minimize( total_energy, x0=x0, bounds = (rmbounds,tbounds),args=(dome_zones[1],phisat, khat, fs, k, 1, 0, 0, 0, cytC_zones[1]),method = 'L-BFGS-B', options = {'ftol':1e-12} )
            theta = ret.x[1]
            rm = ret.x[0] * 1e-6
            print( full_potential(dome_zones[1], shape_rm[1], shape_theta[1], cytC_zones[1], memC_zones[1]) )
            print( "nernst: ", chem_pot_calc(memC_zones[1], cytC_zones[1]))
          print(Fore.WHITE + "")
              
          ax.plot(time_list, np.asarray(globals()["internKf" + str(1) + "_list"]) / Kb, 'b', linestyle = '--', label = 'MOOSE internal')
          if int(num_spines_v.text_disp._text) > 1:
            ax.plot(time_list, np.asarray(globals()["internKf" + str(2) + "_list"]) / Kb, 'r', linestyle = '--', label = 'MOOSE internal')

       print("LENGTH ZONES: ", length_zones)
       print("DOME ZONES: ", dome_zones)
       print("CYT CONCS: ", cytC_zones)
       print("MEM CONCS: ", memC_zones)





       print("Dome volume distribution: ", dome_volume_distro)
       print("MEM ZONES: (C) ", memC_zones)
       allot_conc(num_voxels, marker, memC_zones, dome_starts, dome_ends, 'memConc')
  
       allot_volume(num_voxels, marker, dome_starts, dome_ends, 'memv', spine_length_distro, dome_volume_distro)
           
       
       #theta, Hmean and Sphi are assigned to the voxels of respective pools.
       allot_geom_par(num_voxels, dome_starts, dome_ends, shape_rm, shape_rp, shape_theta, shape_phi, marker)

       output_formatting = [round(dome_zones[s],1) for s in range(int(num_spines_v.text_disp._text))]
       if disable_gui == False:
           out_text.set_val("Spine Content \n " + str(output_formatting) + " \n Sum Spine: " + str(sum(dome_zones)) + "\n SUM MEMBRANE: " + str(sum(membrane)) + "\n Total Molecules: " + str(sum(membrane) + sum(cytosol)) )
       else:
           print(output_formatting)

       #Check mass conservation:
      
       if np.abs(sum(membrane) + sum(cytosol)  - int(strength_v.text_disp._text) ) > 0.01:
           print("Mass conservation error: ")
           break

       #Runs Rdesigneur model for dt. This will update the membrane molecules

       moose.start(dt)
       currentSimTime = moose.element('/clock').currentTime
       run_time =  currentSimTime
       print("Current Time: ", run_time)
       
       

       if r % plot_interval == 0:
            cyt_save_temp = moose.vec("/model/chem/dend/cytosol").n
            total_cyt = [0] * len(cyt_save_temp)
            for c in range(len(cyt_save_temp)):
                total_cyt[c] = cyt_save_temp[c]
            cyt_save.append(total_cyt)

       dome_starts = []
       dome_ends = []
       x = []
       y = []
       x2 = np.asarray(x)
       y2 = np.asarray(y)

       #The geometries and their locations need to be recalculated from the updated phitot. This is done here.
       x2, y2, f_m_val, dome_starts, dome_ends, new_dome_area, phi, shape_rm, shape_rp, shape_theta, shape_phi, prev_xbend, x_positions  = distributing_shape(x0, rmbounds, tbounds, nbounds, dome_zones, int(num_spines_v.text_disp._text), spacing, x2, y2, 1, x_positions_initial, run_time, cytC_zones)

       if r % plot_interval == 0:
           x_save.append(x2)
           y_save.append(y2)
      
                 
       if r % plot_interval == 0:
         time_list.append(currentSimTime)
         Lat_diff_time = (0.5 * dendDia)**2 / (2 * Lat_diff)
         ks = moose.element( '/model/chem/dend/ksolve' )
         rates2 = ks.rateVec['/model/chem/dend/reacF']
         internKfs = list(np.where(np.array(rates2) > 1e-12)[0])
         internKfs.append(0)  #padding so that index dont run over
         tempKf = []
         for ikf in range(len(internKfs) - 1):
             if internKfs[ikf] + 1 != internKfs[ikf + 1]:
                tempKf.append(rates2[internKfs[ikf]])
         for i in range(int(num_spines_v.text_disp._text)):
             globals()["Slength" + str(i+1) + "_list"].append(spine_length_distro[i])
             globals()["area_" + str(i + 1) + "_list"].append( dome_zones[i] / ( phisat * shape_phi[i] ) )
             globals()["phi_" + str(i + 1) + "_list"].append( shape_phi[i] )
             globals()["phitot" + str(i + 1) + "_list"].append( dome_zones[i] )
             globals()["rd" + str(i + 1) + "_list"].append( shape_rm[i] )
             globals()["mu" + str(i + 1) + "_list"].append( f_cm * (1/Lat_diff_time) * np.exp(-chem_pot_partition(dome_zones[i], 1, cytC_zones[i], shape_rm[i] , shape_theta[i])) )
             globals()["internKf" + str(i + 1) + "_list"].append( tempKf[i] )

    if disable_gui == False:
      ax1.plot(x2,y2, color = 'k', label = 'optimised')
      ax1.plot(prev_X,prev_Y, color = 'b', label = 'Initial')
      if max(prev_Y) > max(y2):
         ax1.set_ylim(0,max(prev_Y)+0.05)
      else:   
         ax1.set_ylim(0,max(y2)+0.05)
      ax1.set_xlim(0, Length)
      ax1.set_xlabel("X",fontsize = 18)
      ax1.legend()
    num_seg = len(moose.vec("/model/chem/dend/cytosol").n)
    cytosol = moose.vec("/model/chem/dend/cytosol").n
    membrane = moose.vec("/model/chem/dend/membrane").n
    #Check Length
    x_vec = np.linspace(0,Length,num_seg)  #Check

    time_data = []
    tdf_labels = []
    time_data.append(time_list)
    tdf_labels.append("Time")
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["phi_" + str(i + 1) + "_list"])
        tdf_labels.append("phi_" + str(i + 1))
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["phitot" + str(i + 1) + "_list"])
        tdf_labels.append("phitot" + str(i + 1))
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["rd" + str(i + 1) + "_list"])
        tdf_labels.append("Rd" + str(i + 1))
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["area_" + str(i + 1) + "_list"])
        tdf_labels.append("A_" + str(i + 1))
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["mu" + str(i + 1) + "_list"])
        tdf_labels.append("MU" + str(i + 1))
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["Slength" + str(i + 1) + "_list"])
        tdf_labels.append("Slength" + str(i + 1))
    for i in range(int(num_spines_v.text_disp._text)):
        time_data.append(globals()["internKf" + str(i + 1) + "_list"])
        tdf_labels.append("internKf" + str(i + 1))

    #time_data = [ time_list, phitot_list, phi_1_list, phi_2_list, phi_3_list, area_1_list, area_2_list, area_3_list, phitot1_list, phitot2_list, phitot3_list, rd1_list, rd2_list, rd3_list, mu1_list, mu2_list, mu3_list]
    tdf = pd.DataFrame(time_data)
    tdf = tdf.transpose()
    tdf.columns = tdf_labels
    #tdf.columns = ["Time", "phitot", "phi_1", "phi_2", "phi_3", "A_1", "A_2", "A_3", "phitot1", "phitot2","phitot3", "Rd1", "Rd2", "Rd3", "MU1", "MU2", "MU3"]
    tdf.to_csv("N_" + str(int(num_spines_v.text_disp._text)) + "_timeData_S_" + str(int(strength_v.text_disp._text)) + "_Spacing_" + str(float(space_v.text_disp._text)) + "_Scale_" + str(scales[0]) + ".csv", index = False)

   # if disable_gui == True:
   #   if int(sys.argv[5]) == 0:
   #      tdf.to_csv("N_" + str(int(num_spines_v.text_disp._text)) + "_timeData_S_" + str(int(strength_v.text_disp._text)) + "_Spacing_" + str(float(space_v.text_disp._text)) + "_Scale_" + str(scales[0]) + ".csv", index = False)
    



    if disable_gui == False:
     
       """ 
       ax.set_title("Energy - bending", fontsize=18)
    
       ax.plot(time_list, bending_energy_list)
       ax.set_xlabel("t (s)",fontsize = 18)
       ax.legend()
       """
       ks = moose.element( '/model/chem/dend/ksolve' )
       rates2 = ks.rateVec['/model/chem/dend/reacF']
       #ax.plot(time_list, globals()["mu" + str(1) + "_list"], 'b', label = 'Calculation')
       ax.set_xlabel("t (s)",fontsize = 18)
       ax.legend()
       
       ax3.clear()
       ax3.set_title("Mem Vs Cyt")
       ax3.set_ylabel("Cytosol",fontsize=18,color = 'b')
       ax3.set_xlabel("X", fontsize = 18)

    x_vec = np.linspace(0,Length,num_seg)  #Check
    
    if disable_gui == False:
      ax4.clear()
      ax4.set_ylabel("Membrane",fontsize=18, color = 'r')
      ax4.plot(x_vec,moose.vec("/model/chem/dend/membrane").n,color = 'r')
      ax5.plot(x_vec,np.asarray(rates2) / Kb, color = 'r')
      ax5.set_title("$k_{f}$")
      ax5.set_xlabel("x $\mu m$")

    if disable_gui == False:    
      ax4.set_ylim(ymin = 0)
      ax3.plot(x_vec, moose.vec("/model/chem/dend/cytosol").n, color = 'b')
    print("WRITING XML")
    writeXML_td(cyt_save,"Cyt_N_" + str(int(num_spines_v.text_disp._text)) + "_timeData_S_" + str(int(strength_v.text_disp._text)) + "_Spacing_" + str(float(space_v.text_disp._text)) + "_Scale_" + str(scales[0]) + ".xml")   
    writeXML_td(x_save,"X_" + str(int(num_spines_v.text_disp._text)) + "_timeData_S_" + str(int(strength_v.text_disp._text)) + "_Spacing_" + str(float(space_v.text_disp._text)) + "_Scale_" + str(scales[0]) + ".xml")   
    writeXML_td(y_save,"Y_" + str(int(num_spines_v.text_disp._text)) + "_timeData_S_" + str(int(strength_v.text_disp._text)) + "_Spacing_" + str(float(space_v.text_disp._text)) + "_Scale_" + str(scales[0]) + ".xml")   


def dendShape( theta, rp, rm,  n, prev_Lflat, x_positions_initial, run_time):
    dth = theta / 100
    xbend = 2 * np.sin( theta ) * (rp + rm )

    if run_time == 0:
      if n == 0:
         Lflat = prev_Lflat - xbend/2.0
      else:
         Lflat = prev_Lflat
    else:
      Lflat = x_positions_initial[n] - xbend/2.0

    x = [0.0, Lflat]
    y = [0.0, 0.0]
    L_truncate = 2 * Lflat + xbend
    
    # To show centre of curvature

    """
    x.extend( [Lflat, Lflat] )
    y.extend( [rp, 0.0] )

    """
    #Starting from flat through the saddle
    x.extend( [ Lflat + rp * np.sin( th + dth ) for th in np.arange(0.0, theta * 0.999999, dth ) ] )
    y.extend( [ rp * (1-np.cos( th + dth) ) for th in np.arange(0.0, theta * 0.999999, dth ) ] )
    dome_start = x[-1]*1e-6
    #End of the saddle
    xlast = x[-1]
    ylast = y[-1]
    #Going to the middle of dome
    xoffset = rm * np.sin( theta ) + xlast
    yoffset = -rm * np.cos( theta ) + ylast

    #Going from the beginning of the dome to the middle.
    x.extend( [ xoffset - rm * np.sin( th ) for th in np.arange (theta, 0, -dth ) ] )
    y.extend ([yoffset + rm * np.cos( th ) for th in np.arange (theta, 0, -dth ) ] )
    xlast = x[-1]
    ylast = y[-1]
    x.extend( [ L_truncate - i for i in x[::-1] ] )
    y.extend( y[::-1] )
    dome_end_x = xoffset + rm * np.sin( theta )
    dome_end = x[-1]
    xlast = x[-1]
    ylast = y[-1]

    x_pos = y.index(max(y))
    return np.array( x ), np.array( y ), dome_start, dome_end_x * 1e-6, Lflat, x[x_pos]


if disable_gui == False:
  fig = plt.figure(figsize=(18,20))
  plt.subplots_adjust(bottom=0.05, right=0.76, top=0.92, hspace = 1.0, wspace = 0.45)
  ax = fig.add_subplot(5, 2, 1)

global l
global l1
global l2
global l3
global X
global Y
global Z
global ax_dome
global dome_v
global dome_v2
global steady_mem
global steady_cyto
marking.action = True

if disable_gui == False:
  ax.set_xlabel("$\Theta$")
  ax.set_ylabel("$Energy$")
  ax1 = fig.add_subplot(5, 2, 2)
  ax2 = fig.add_subplot(5, 2, 3)
  ax2_1 = ax2.twinx()
  ax3 = fig.add_subplot(5, 2, 4)
  ax4 = ax3.twinx()
  ax5 = fig.add_subplot(5, 2, 5)
  ax6 = fig.add_subplot(5, 2, 6)
  ax9 = fig.add_subplot(5, 2, 7)
  ax10 = fig.add_subplot(5, 2, 8)
  ax11 = fig.add_subplot(5, 2, 9)
  ax12 = fig.add_subplot(5, 2, 10)

  ax_num_spines = plt.axes([0.93, 0.1, 0.03, 0.03])
  ax_theta2 = plt.axes([0.93, 0.18, 0.03, 0.03])
  ax_runtime = plt.axes([0.93, 0.25, 0.03, 0.03])
  ax_start = plt.axes([0.78, 0.45, 0.2, 0.05])
  ax_text = plt.axes([0.79, 0.85, 0.2, 0.1])
  ax_dome = plt.axes([0.93, 0.4, 0.03, 0.03])
  ax_dome2 = plt.axes([0.93, 0.35, 0.03, 0.03])
  ax_spacing = plt.axes([0.93, 0.3, 0.03, 0.03])
  ax_strength = plt.axes([0.93, 0.05, 0.03, 0.03])

  dome_v = TextBox(ax_dome, 'Scaling Dome 1  ', initial = ".5 0.0")
  dome_v.label.set_size(18)
  dome_v2 = TextBox(ax_dome2, 'Scaling Dome 2  ', initial = "1")
  dome_v2.label.set_size(18)
  num_spines_v = TextBox(ax_num_spines, 'num spines  ', initial = "3")
  num_spines_v.label.set_size(18)
  theta_v2 = TextBox(ax_theta2, 'Unused  ', initial = "0.5")
  theta_v2.label.set_size(18)
  runtime_v = TextBox(ax_runtime, 'runtime (s)  ', initial = "1")
  runtime_v.label.set_size(18)
  start_b = Button(ax_start, "What's the answer?",color = "green")
  start_b.label.set_size(18)
  out_text = TextBox(ax_text, "")
  out_text.label.set_size(18)

  space_v = TextBox(ax_spacing, 'spacing $\mu m$  ', initial = "6")
  space_v.label.set_size(18)
  strength_v = TextBox(ax_strength, 'strength', initial = "8000")
  strength_v.label.set_size(18)
  scaling_factor = 0.2

  scale_string = str(scaling_factor)    
  print(int(num_spines_v.text_disp._text))
  for i in range(1, int(num_spines_v.text_disp._text)):
      scale_string = scale_string + " "
      print((1 - scaling_factor) / ( int(num_spines_v.text_disp._text) - 1 ))
      scale_string = scale_string + str( (1 - scaling_factor) / ( int(num_spines_v.text_disp._text) - 1 ) )
  print(scale_string)            

  if int(num_spines_v.text_disp._text) == 3:
    scales_distribution = [scaling_factor, 1 - 2 * scaling_factor  , scaling_factor]   
    scale_string = str(scales_distribution[0])
    for i in range(1, int(num_spines_v.text_disp._text)):
      scale_string =  scale_string + " " + str(scales_distribution[i])
    print("SCALE STRING: ", scale_string)  
  dome_v.set_val(scale_string)
  
  start_b.on_clicked(update)
  plt.show()
else:
  ax_num_spines = plt.axes([0.93, 0.1, 0.03, 0.03])
  num_spines_v = TextBox(ax_num_spines, 'num spines  ', initial = "5")
  num_spines_v.set_val(str(sys.argv[4]))
  ax_strength = plt.axes([0.93, 0.05, 0.03, 0.03])
  strength_v = TextBox(ax_strength, 'strength', initial = "40000")
  strength_v.set_val(str(sys.argv[6]))
  ax_dome = plt.axes([0.93, 0.4, 0.03, 0.03])
  dome_v = TextBox(ax_dome, 'Scaling Dome 1  ', initial = "0.65 0.35")
  scaling_factor = float(sys.argv[1])
  if int(num_spines_v.text_disp._text) == 1:
      dome_v.set_val(str(scaling_factor) + " " + str(0))
  scale_string = str(scaling_factor)    
  print(int(num_spines_v.text_disp._text))
  for i in range(1, int(num_spines_v.text_disp._text)):
      scale_string = scale_string + " "
      print((1 - scaling_factor) / ( int(num_spines_v.text_disp._text) - 1 ))
      scale_string = scale_string + str( (1 - scaling_factor) / ( int(num_spines_v.text_disp._text) - 1 ) )
  print(scale_string)            
  
  if int(num_spines_v.text_disp._text) == 3:
    scales_distribution = [scaling_factor, 1 - 2 * scaling_factor  , scaling_factor]   
    scale_string = str(scales_distribution[0])
    for i in range(1, int(num_spines_v.text_disp._text)):
      scale_string =  scale_string + " " + str(scales_distribution[i])
    print("SCALE STRING: ", scale_string)  
  dome_v.set_val(scale_string)

  ax_spacing = plt.axes([0.93, 0.3, 0.03, 0.03])
  space_v = TextBox(ax_spacing, 'spacing $\mu m$  ', initial = "2")
  space_v.set_val("2")
  ax_runtime = plt.axes([0.93, 0.25, 0.03, 0.03])
  runtime_v = TextBox(ax_runtime, 'runtime (s)  ', initial = "1")
  runtime_v.set_val("1")
  update(1)




    



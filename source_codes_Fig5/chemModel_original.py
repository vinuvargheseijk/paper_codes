import totalEnergy
import moose
import numpy as np
import rdesigneur as rd
import tag
import xml.etree.ElementTree as ET
import pandas as pd
#import matplotlib.pyplot as plt

k = totalEnergy.k
khat = totalEnergy.khat
Na = totalEnergy.Na
RT = totalEnergy.RT
mu0 = totalEnergy.calc_mu0()
KBT = totalEnergy.KBT
Cp = totalEnergy.Cp
diffL = totalEnergy.diffL
phisat = totalEnergy.phisat
k_agg = totalEnergy.k_agg
dendDia = totalEnergy.dendDia
Length = totalEnergy.Length
part_coeff = totalEnergy.part_coeff
Kb = tag.Kb
diffDt = 0.0001
chemDt = diffDt
cyt_conc = tag.cyt_conc
comptLen = tag.comptLen * 1e-6
elec_dx = comptLen
num_segments = int( Length * 1e-6 / elec_dx )
#cyt_motor_rate = tag.cyt_motor_rate #m/s


def dummy_stimulus_delivery():
  ECD = 0.3e-12
  stim1=moose.vec( '/model/chem/dend/stim1' ).n
  comparts = moose.wildcardFind( '/model/elec/dend#' )
  gri = 0
  runtime = tag.runtime
  tabDt = chemDt
  size = int( runtime / tabDt )
  for compt in comparts:
     tab = moose.element( compt.path + "/glu/tab" )
     moose.connect( tab, "output", compt.path + "/NMDA", "activation" )
     tab.startTime = 0
     tab.stopTime = runtime
     tab.stepSize = 0  #Uses the current time
     tab.vector = [1e-3] * size

def stimulus_delivery():
  df = pd.DataFrame()  
  #ECD = 0.09e-12
  ECD = 0.3e-12
  stim_type = tag.stim_type
  stim1=moose.vec( '/model/chem/dend/stim1' ).n
  stim_times = tag.glut_stim_times
  stim_locations = tag.glut_stim_locations
  stim_y = tag.glut_stim_y
  comparts = moose.wildcardFind( '/model/elec/dend#' )
  gri = 0
  runtime = tag.runtime
  tabDt = chemDt
  size = int( runtime / tabDt )
  trunc_size = int( stim_times[-1] / tabDt ) + int( 2 / tabDt )
  total_molecules = 0
  dt_b = 0.025  #Burst time step
  sflag1 = []
  lflag1 = []
  yflag1 = []
  burst_num_pulse = tag.burst_num_pulse
  for scount in range( len( stim_times ) ):
      if stim_type[scount]["stim"] == 'b':
         print("RECORDING BURST") 
         for i in range( burst_num_pulse ):
             sflag1.append( stim_times[scount] + i * dt_b )     
             lflag1.append( stim_locations[scount] )
             yflag1.append( stim_y[scount] )
      if stim_type[scount]["stim"] == 's':
             sflag1.append( stim_times[scount] ) 
             lflag1.append( stim_locations[scount] ) 
             yflag1.append( stim_y[scount] )
  print("STIM TIMES in CHEM: ", sflag1)           
  if tag.bg == True:
      for bcount in range( int(tag.num_bg_times) ):
        if bg_stim_type == 'b':
            for i in range( burst_num_pulse ):
                sflag1.append( tag.bg_start + bcount * tag.bg_dt + i * dt_b )     
                lflag1.append( tag.bg_locs[bcount] )
                yflag1.append( y_pos )
        if bg_stim_type == 's':
                sflag1.append( tag.bg_start + bcount * tag.bg_dt )
                lflag1.append( tag.bg_locs[bcount] )
                yflag1.append( y_pos )
  for compt in comparts:
     rec_temp = []
     tab = moose.element( compt.path + "/glu/tab" )
     moose.connect( tab, "output", compt.path + "/NMDA", "activation" )
     tab.startTime = 0
     tab.stopTime = runtime
     tab.stepSize = 0  #Uses the current time
     for i in range( size + 1000 ): 
        c_time = i * tabDt
        if c_time < sflag1[-1] + 2:
          for j in range( len( sflag1 ) ):
            #Check
            #if c_time - sflag1[j] <= 0.0:
            #    stim1[gri] = 0.0   #The stimulus at t = 0  is zero
            if c_time - sflag1[j] > 0:
                if c_time <  tag.end_stim_time and tag.bg == True:
                       Source = tag.bg_int
                if tag.bg == False:   
                       Source = tag.Source
                conc_factor = 1 / (totalEnergy.Na * yflag1[j])
                deltaZ = yflag1[j]
                deltaY = elec_dx
                #stim1[gri] += conc_factor * (Source/(4*np.pi*(c_time-sflag1[j])*1e3*ECD))*np.exp(- (gri * elec_dx - lflag1[j])**2 /(4*ECD*(c_time-sflag1[j])*1e3))*np.exp(-(yflag1[j]**2)/(4*ECD*(c_time-sflag1[j])*1e3))
                stim1[gri] += (elec_dx * deltaY * totalEnergy.Na * deltaZ)  * conc_factor * (Source/(4*np.pi*(c_time-sflag1[j])*1e3*ECD))*np.exp(- (gri * elec_dx - lflag1[j])**2 /(4*ECD*(c_time-sflag1[j])*1e3))*np.exp(-(yflag1[j]**2)/(4*ECD*(c_time-sflag1[j])*1e3))
                if stim1[gri] > 1e3:
                    print("STIM EXCESS AT: ", c_time, sflag1[j])
        else:
            stim1[gri] = 0.0
        rec_temp.append( stim1[gri] )
        stim1[gri] = 0
     #df[str(gri)] = rec_temp   
     gri += 1
     tab.vector = rec_temp
     print("LFLAG and SFLAG: ", lflag1, sflag1)          
     print("YFLAG: ", yflag1)          
     #plt.plot(tab.vector)
     #plt.ylim(-10,100)
     total_molecules = total_molecules + sum(tab.vector)
     #print(tab.vector)   

  #min_imax = 1e6
  #mid_p = int(len(comparts) / 2.0)
  #for sl in [mid_p -1, mid_p, mid_p + 1]:
  #    imax = df[str(sl)].argmax()
  #    if imax < min_imax:
  #        min_imax = imax
  #sum_g = 0
  #for i in range(len(comparts)):
  #     sum_g = sum_g + df[str(i)][min_imax]
  #print("Sum of molecules: ", sum_g)     
  #plt.show()

def chem_model():
    moose.Neutral( '/library' )
    meshName = '/library/cm'
    moose.loadModel('./acc49_IR_v4_high_RT.g', meshName, 'ee')
    #moose.SBML.mooseReadSBML('./acc49_v22.xml', meshName)

    #moose.CubeMesh( '/library/cm/dend' )
    Ca=moose.element(meshName+'/kinetics/Ca')
    Rac_GTP=moose.element(meshName+'/kinetics/Rac_GTP')
    Rac_GDP=moose.element(meshName+'/kinetics/Rac_GDP')
    Tiam_Rev_Reac = moose.element(meshName+'/kinetics/Tiam_Rev_Reac')
    #Enz_tot_CaM_CaMKII = moose.element(meshName+'/kinetics/CaMKII_gr/CaMKII_thr286/Enz_tot_CaM_CaMKII')
    tot_act_CaMKII = moose.element(meshName+'/kinetics/CaMKII_gr/tot_act_CaMKII')
    tot_CaM_CaMKII = moose.element(meshName+'/kinetics/CaMKII_gr/tot_CaM_CaMKII')
    tot_autonomous_CaMKII = moose.element(meshName+'/kinetics/CaMKII_gr/tot_autonomous_CaMKII')
    Enz_Tiam = moose.element(meshName+'/kinetics/Tiam1/Enz_Tiam')
    Tiam_Rev_Reac = moose.element(meshName+'/kinetics/Tiam_Rev_Reac')
    cyt_act_reac = moose.element(meshName+'/kinetics/cyt_act_reac')
    Rac_Rev_Reac = moose.element(meshName+'/kinetics/Rac_Rev_Reac')
    Tiam_Metab_Reac = moose.element(meshName+'/kinetics/Tiam_Metab_Reac')
    Tiam_Metab = moose.element(meshName+'/kinetics/Tiam_Metab')
    Rac_Metab_Reac = moose.element(meshName+'/kinetics/Rac_Metab_Reac')
    Rac_Metab = moose.element(meshName+'/kinetics/Rac_Metab')
    marker = moose.Pool( meshName + '/kinetics/marker' )
    pot_eval = moose.Pool( meshName + '/kinetics/pot_eval' )
    Hmean = moose.Pool( meshName + '/kinetics/Hmean' )
    Sphi = moose.Pool( meshName + '/kinetics/Sphi' )
    cytosol = moose.element( meshName + '/kinetics/cytosol' )
    dimer = moose.element( meshName + '/kinetics/IRSp53_dimer' )
    memConc = moose.Pool( meshName + '/kinetics/memConc' )
    cytConc = moose.Pool( meshName + '/kinetics/cytConc' )
    memv = moose.Pool( meshName + '/kinetics/memv' )
    membrane = moose.Pool( meshName + '/kinetics/membrane' )
    theta = moose.Pool( meshName + '/kinetics/theta' )
    pot_nernst = moose.Pool( meshName + '/kinetics/pot_nernst' )
    saved_kf = moose.Pool( meshName + '/kinetics/saved_kf' )
    stim1 = moose.Pool( meshName + '/kinetics/stim1' )
    non_diff_mols = [marker, membrane, pot_eval, Hmean, Sphi, memConc, cytConc, memv, theta, pot_nernst, saved_kf, stim1, Tiam_Metab, Rac_Metab, tot_act_CaMKII, tot_CaM_CaMKII, tot_autonomous_CaMKII]
    print("Printing concs:")
    print("Tiam_inact: ", moose.element(meshName+'/kinetics/Tiam_inact').conc)
    print("Rac_GDP: ", moose.element(meshName+'/kinetics/Rac_GDP').conc)

    reacF = moose.Reac( meshName + '/kinetics/reacF' )
    reacB = moose.Reac( meshName + '/kinetics/reacB' )


    moose.connect( reacF, 'sub', cytosol, 'reac')
    moose.connect( reacF, 'prd', membrane, 'reac')
    moose.connect( reacB, 'sub', membrane, 'reac')
    moose.connect( reacB, 'prd', cytosol, 'reac')

    reacF.Kf = 0
    reacF.Kb = 0
    reacB.Kf = 0
    reacB.Kb = 0
    func_kf = moose.Function(meshName + '/kinetics/func_kf')
    func_pot_eval = moose.Function(meshName + '/kinetics/func_pot_eval')
    func_pot_nernst_eval = moose.Function(meshName + '/kinetics/func_pot_nernst_eval')
    saved_kf_func = moose.Function( meshName + '/kinetics/saved_kf_func' )


    func_pot_eval.x.num = 4
    func_pot_nernst_eval.x.num = 3
    func_kf.x.num = 3
    saved_kf_func.x.num = 3

    voxel_volume = np.pi * (0.5 * dendDia)**2 * diffL
    func_kf.expr = str(part_coeff * Kb) + " * x0 * exp( -1 * (x1 / " + str(RT) + ") ) + 0 * x2"   #Note: mu0 and RT are in J/molecule

    print(func_kf.expr)

    saved_kf_func.expr = str(part_coeff * Kb) + " * x0 * exp( -1 * (x1 / " + str(RT) + ") ) + 0 * x2"   #Note: mu0 and RT are in J/molecule


    func_pot_eval.expr = "x0 * (" + str( 0.5 * khat / phisat ) + " * (2 * x1 - " + str(Cp) + ")^2 + " + str(KBT) + " * ( log(x2 + (x2 < 1e-2) * 1) - log(1 - x2) ) + " + str(KBT) + " * ( 1 + " + str(totalEnergy.k_agg) + " * x2 + " + str(totalEnergy.k_agg2) + " * 3 * x2^2 ) + x3 )" 

    func_pot_nernst_eval.expr = "x0 * ( " + str(mu0) + " + " + str(RT) + " * log( (x1 + 1e-32) / (x2 + 1e-32) ) )" #This should give unit in Joules / molecule

    moose.connect( marker, 'nOut', func_kf.x[0], 'input' ) 
    moose.connect( pot_eval, 'nOut', func_kf.x[1], 'input' )
    moose.connect( Sphi, 'nOut', func_kf.x[2], 'input' )
    moose.connect( func_kf, 'valueOut', reacF, 'setNumKf' )
    
    moose.connect( marker, 'nOut', saved_kf_func.x[0], 'input' ) 
    moose.connect( pot_eval, 'nOut', saved_kf_func.x[1], 'input' )
    moose.connect( Sphi, 'nOut', saved_kf_func.x[2], 'input' )
    moose.connect( saved_kf_func, 'valueOut', saved_kf, 'setN' )
    
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

    if tag.chemTurnOff:
       cyt_act_reac.Kf = 0
       cyt_act_reac.Kb = 0
    else:
       Tiam_Metab_Reac.Kf = tag.Tiam_Metab_Reac_Kf
       Tiam_Metab_Reac.Kb = tag.Tiam_Metab_Reac_Kf
       Rac_Metab_Reac.Kf = tag.Rac_Metab_Reac_Kf
       Rac_Metab_Reac.Kb = tag.Rac_Metab_Reac_Kf
       #Enz_tot_CaM_CaMKII.Km = tag.Enz_tot_CaM_CaMKII_Km
       #Enz_Tiam.Km = tag.Enz_Tiam_Km
       #Enz_Tiam.kcat = tag.Enz_Tiam_kcat
       #Enz_tot_CaM_CaMKII.kcat = tag.Enz_tot_CaM_CaMKII_kcat

       Rac_Rev_Reac.Kf = tag.Rac_Rev_Reac_Kf
       Rac_Rev_Reac.Kb = 0.00001

       Tiam_Rev_Reac.Kf = tag.Tiam_Rev_Reac_Kf
       Tiam_Rev_Reac.Kb = 0.00001



    file_molWt='prot_wt.xml'
    den=6.0*3.14*np.cbrt(3.0/(4.0*3.14*0.73*6.02*1e23))
    num=1e3*1e6*1.38*1e-23
    co=num/den
    df_diff = pd.DataFrame()
    molNames = []
    molWts = []
    molDiffs = [] 
    print('diffConst updated')
    for prot in moose.wildcardFind(meshName+'/kinetics/#[ISA=Pool]'):
        if prot not in non_diff_mols: 
            mol_name=ET.parse(file_molWt).find(str(prot.name))
            print('(Kinetics: The name is ')
            print(prot)
            print(mol_name)
            molWt=mol_name.text.split()
            print(molWt[0])
            molNames.append(molWt[0])
            molWts.append(molWt[1])
            invmol=1/np.cbrt(float(molWt[1])*1e3)
            moose.element(prot).diffConst=(21.83482335*invmol-0.14591424)*1e-12
            print(moose.element(prot).diffConst)
            molDiffs.append(moose.element(prot).diffConst * 1e12)
    for prot in moose.wildcardFind(meshName+'/kinetics/CaMKII_gr/#[ISA=Pool]'):
        if prot not in non_diff_mols: 
            mol_name=ET.parse(file_molWt).find(str(prot.name))
            print('CaMKII: The name is ')
            print(prot)
            print(mol_name)
            molWt=mol_name.text.split()
            invmol=1/np.cbrt(float(molWt[1])*1e3)
            moose.element(prot).diffConst=(21.83482335*invmol-0.14591424)*1e-12
            print(moose.element(prot).diffConst)
            molNames.append(molWt[0])
            molWts.append(molWt[1])
            molDiffs.append(moose.element(prot).diffConst * 1e12)
    df_diff["name"] = molNames
    df_diff["wt"] = molWts
    df_diff["D"] = molDiffs
    df_diff.to_csv("./DiffC.csv")
    Ca.diffConst=20e-12
    print(Ca.diffConst)

    print('DIFFUSION FOR IRSp53 GROUP')
 
    invmol = 1/np.cbrt(22.0*1e3)

    Rac_GTP.diffConst = (21.83482335*invmol-0.14591424)*1e-12
    print("RAC diff: ", Rac_GTP.diffConst)
    Rac_GDP.diffConst = Rac_GTP.diffConst

    membrane.diffConst = tag.mem_diffConst   #Check
    cytosol.concInit = cyt_conc
    dimer.concInit = tag.dimer_conc

def rdes_build():
    cm = chem_model()
    rdes = rd.rdesigneur(
    useGssa = False,        
    turnOffElec = False,
    combineSegments = False,
    stealCellFromLibrary = True,
    diffusionLength = diffL,
    diffDt = diffDt,
    chemDt = chemDt,
    #cellProto = [['soma', 'soma']],  #Check
    cellProto = [['ballAndStick', 'soma', 10e-6, 10e-6, dendDia, (Length) * 1e-6, num_segments]],  #Check
    chemProto = [['cm', 'cm']],
    chemDistrib = [['cm', 'dend#', 'install', '1' ]],    
    chanProto = [['make_Na()', 'Na'],
                 ['make_K_DR()', 'K_DR'],
                #['make_K_A()', 'K_A' ],
                 ['make_NMDA()','NMDA'],
                 ['make_Ca_conc()', 'Ca_conc'],
                 ['make_glu()', 'glu' ],
                 ['make_Ca()', 'Ca' ],
                ],
    chanDistrib = [['Na', 'soma', 'Gbar', '30' ],
                   ['K_DR', 'soma', 'Gbar', '25' ],
                   ['NMDA', 'dend#', 'Gbar', '9' ],  #Check
                   ['Ca_conc', 'dend#', 'tau', str(tag.CaTau)], #Check
                   ['glu', 'dend#', 'Gbar', '3' ], #Check
                   ['Ca', 'dend#', 'Gbar', '0.0' ],
                 ],
    passiveDistrib = [['soma', 'CM', '0.01', 'Em', '-0.06']],
    adaptorList = [
                   [ 'Ca_conc', 'Ca', 'dend/Ca', 'conc', 0.00008, tag.CaScale ]
                  ],
    )
    gtab = moose.StimulusTable( '/library/glu/tab' )
    moose.connect( gtab, 'output', '/library/glu', 'activation' )
    print("Building model")
    model = rdes.buildModel('/model')
    print(moose.le('/model/elec'))
    #moose.setClock(30, diffDt)   #PyRun clock tick 30
    #moose.setClock(10, chemDt)
    #for i in range( 11, 18 ):
    #   moose.setClock( i, chemDt )
    #moose.setClock( 18, chemDt )


    
rdes_build()
print("After build: ", moose.le('/model/chem/dend'))
print("Length: ", len(moose.vec('/model/chem/dend/cytosol').n))
if tag.chemTurnOff == False:
   stimulus_delivery()
else:
   dummy_stimulus_delivery()
moose.reinit()


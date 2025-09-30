import totalEnergy
import moose
import numpy as np
import rdesigneur as rd
import tag

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
diffDt = 0.01
cyt_conc = 1e-5
#cyt_motor_rate = tag.cyt_motor_rate #m/s

def chem_model():
    moose.Neutral( '/library' )
    moose.Neutral( '/library/cm' )
    moose.CubeMesh( '/library/cm/dend' )
    marker = moose.Pool( '/library/cm/dend/marker' )
    pot_eval = moose.Pool( '/library/cm/dend/pot_eval' )
    Hmean = moose.Pool( '/library/cm/dend/Hmean' )
    Sphi = moose.Pool( '/library/cm/dend/Sphi' )
    cytosol = moose.Pool( '/library/cm/dend/cytosol' )
    memConc = moose.Pool( '/library/cm/dend/memConc' )
    cytConc = moose.Pool( '/library/cm/dend/cytConc' )
    memv = moose.Pool( '/library/cm/dend/memv' )
    membrane = moose.Pool( '/library/cm/dend/membrane' )
    theta = moose.Pool( '/library/cm/dend/theta' )
    pot_nernst = moose.Pool( '/library/cm/dend/pot_nernst' )
    saved_kf = moose.Pool( '/library/cm/dend/saved_kf' )


    reacF = moose.Reac( '/library/cm/dend/reacF' )
    reacB = moose.Reac( '/library/cm/dend/reacB' )

    cytosol.diffConst = tag.cyt_diffConst
    #cytosol.motorConst = cyt_motor_rate
    membrane.diffConst = tag.mem_diffConst   #Check
    cytosol.concInit = cyt_conc

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
    saved_kf_func = moose.Function( '/library/cm/dend/saved_kf_func' )
    
    


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

def makeTaper():
        comptLen = tag.comptLen * 1e-6
        RM = 1.0
        RA = 10.0
        CM = 0.01
        comptDia = tag.bigDia * 1e-6
        taper = tag.taper
        soma = moose.Neuron( '/library/soma' )
        prev = rd.buildCompt( soma, 'soma', RM = RM, RA = RA, CM = CM, dia = comptDia, x=0, dx=comptLen)
        theta = 0
        x = comptLen
        y = 0.0
        numSegments = int(Length * 1e-6 / comptLen)


        for i in range( numSegments - 1 ):
            #dx = comptLen * np.cos( theta )
            #dy = comptLen * np.sin( theta )
            dx = comptLen
            dy = 0
            print(comptDia)
            print(comptDia - i * taper)
            r = np.sqrt( x * x + y * y )
            theta += comptLen / r
            compt = rd.buildCompt( soma, 'soma' + str(i), RM = RM, RA = RA, CM = CM, x = x, y = y, dx = dx, dy = dy, dia = comptDia - i * taper * 1e-6  )
            moose.connect( prev, 'axial', compt, 'raxial' )
            prev = compt
            x += dx
            y += dy
        print("Length constructed: ", x)    
        print("Final compartment dia: ", comptDia - i * taper)

        return soma

moose.Neutral('/library')
makeTaper()    


def rdes_build():
    cm = chem_model()
    rdes = rd.rdesigneur(
    turnOffElec = True,
    diffusionLength = diffL,
    diffDt = diffDt,
    chemDt = diffDt,
    cellProto = [['soma', 'soma']],  #Check
    chemProto = [['cm', 'cm']],
    chemDistrib = [['cm', 'soma#', 'install', '1' ]],    
    )
    print("Building model")
    model = rdes.buildModel('/model')
    print(moose.le('/model/elec'))
    moose.setClock(30, diffDt)   #PyRun clock tick 30
    moose.setClock(10, diffDt)
    for i in range( 11, 18 ):
       moose.setClock( i, diffDt )
    moose.setClock( 18, diffDt )


    
rdes_build()
print("After build: ", moose.le('/model/chem/dend'))
print("Length: ", len(moose.vec('/model/chem/dend/cytosol').n))
moose.reinit()


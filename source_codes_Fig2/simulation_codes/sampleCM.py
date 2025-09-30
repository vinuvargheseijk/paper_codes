import totalEnergy
import moose
import numpy as np
import rdesigneur as rd

k = totalEnergy.k
khat = totalEnergy.khat
RT = totalEnergy.RT
mu0 = totalEnergy.mu0
KBT = totalEnergy.KBT
Cp = totalEnergy.Cp
diffL = totalEnergy.diffL
phisat = totalEnergy.phisat
k_agg = totalEnergy.k_agg
dendDia = totalEnergy.dendDia
Length = totalEnergy.Length
Lat_diff = 1.0e-12
Lat_diff_time = (0.5 * dendDia)**2 / (2 * Lat_diff)
Kb = 1/Lat_diff_time
diffDt = 0.0005

def chem_model():
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

    reacF.Kf = 1
    reacF.Kb = 0
    reacB.Kf = 0
    reacB.Kb = 0

    reacB.Kf = 0.0

    reacB.Kb = 0.0

    reacF.Kb = 0.0

    cytosol.concInit = 0.001



def rdes_build():
    cm = chem_model()
    rdes = rd.rdesigneur(
    turnOffElec = True,
    diffusionLength = diffL,
    diffDt = diffDt,
    chemDt = diffDt,
    cellProto = [['somaProto', 'soma', dendDia, (Length) * 1e-6]],  #Check
    chemProto = [['cm', 'cm']],
    chemDistrib = [['cm', 'soma', 'install', '1' ]],    
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

moose.reinit()
moose.start(1)

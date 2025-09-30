import numpy as np
import tag
import totalEnergy as TE

phisat = tag.phisat
Kb = tag.Kb
def wave_speed_calc(rd, theta, phi, conc, kf):
    area = 2 * np.pi * rd**2 * (1 - np.cos(theta))
    phitot = area * phi * phisat
    L = rd * np.sin( theta )
    num_voxels = 2 * L / TE.diffL #Have to take the entire dome length, L = halfwidth
    CI = np.sqrt( 1 / ( 2 * np.pi * phi * phisat ) )
    speed = 0.5 * phitot**-0.5 * CI * np.sqrt(1 + np.cos(theta)) + CI**2 * np.cos(theta) / L
    nc = conc * np.pi * (tag.dendDia**2 / 4) * TE.Na * TE.diffL
    nm = phitot / num_voxels
    kf_nc = kf * nc
    kb_nm = Kb * nm
    print("LENGTH OF DOME: ", L)
    print("kf_c, kb_m: ", kf_nc, kb_nm)
    return speed * (kf_nc - kb_nm)

import numpy as np
import matplotlib.pyplot as plt
import totalEnergy as TE
import chemModel as CM
import moose


num_voxels = len(moose.vec('/model/chem/dend/cytosol').n)
mem_coords = []
for nv in range(num_voxels):
    mem_coords.append( moose.element('/model/chem/dend/mesh['+str(nv)+']').Coordinates[0] )

def oneD_diffusion_profile(location_um, strength, t, diff, mem_coords):
   strength_factor = ( strength * TE.diffL ) / ( np.sqrt(4 * np.pi * diff * t) )
   location_coord = int( location_um / TE.diffL ) 
   print("STIM LOC: ", location_coord)
   stim_vec = []
   cyt_temp = moose.vec('/model/chem/dend/cytosol').n
   for mc in range(len(mem_coords)):
       diff_factor = np.exp( - (location_um - mem_coords[mc])**2 / (4 * diff * t) )
       stim_vec.append( strength_factor * diff_factor )
   cyt_update = np.asarray(cyt_temp) + np.asarray(stim_vec)    
   moose.vec('/model/chem/dend/cytosol').n = cyt_update
   return stim_vec

def twoD_diffusion_profile(location_um, strength, t, diff, distance, mem_coords):
   strength_factor = ( strength * (TE.diffL)**2 ) / ( 4 * np.pi * diff * t )
   location_coord = int( location_um / TE.diffL )
   print("STIM LOC: ", location_coord)
   stim_vec = []
   cyt_temp = moose.vec('/model/chem/dend/cytosol').n
   for mc in range(len(mem_coords)):
       diff_factor = np.exp( - (location_um - mem_coords[mc])**2 / (4 * diff * t)  - distance**2 / (4 * diff * t) )
       stim_vec.append( strength_factor * diff_factor )
   cyt_update = np.asarray(cyt_temp) + np.asarray(stim_vec)
   moose.vec('/model/chem/dend/cytosol').n = cyt_update
   return stim_vec   

stim_loc = (TE.Length * 1e-6) / 2.0 
print("STIM LOC: ", stim_loc)
times = [1, 5, 10, 20, 30, 50, 70, 90, 100]
for t in times: 
  stim_vec = twoD_diffusion_profile(stim_loc, 8000, t, 0.1e-12, 1e-6, mem_coords)
  print("SUM: ", sum(stim_vec))
  plt.plot(mem_coords, stim_vec)
plt.show()


import numpy as np
import matplotlib.pyplot as plt
import totalEnergy as TE
import chemModel as CM
import moose
import csv

def oneD_diffusion_profile(location_um, strength, t, diff, mem_coords, current_time):
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
   with open("./stim_q.csv", 'a') as stim_file:
       writer = csv.writer(stim_file)
       writer.writerow( [ current_time, sum(stim_vec)] )
       stim_file.close()


def twoD_diffusion_profile(location_um, strength, t, diff, distance, mem_coords):
   strength_factor = ( strength * (TE.diffL)**2 ) / ( 4 * np.pi * diff * t )
   location_coord = int( location_um / TE.diffL ) 
   print("STIM LOC: ", location_coord)
   stim_vec = []
   cyt_temp = moose.vec('/model/chem/dend/cytosol').n
   for mc in range(len(mem_coords)):
       diff_factor = np.exp( - (location_um - mem_coords[mc])**2 / (4 * diff * t)  - (distance - mem_coords[mc])**2 / (4 * diff * t) )
       stim_vec.append( strength_factor * diff_factor )
   cyt_update = np.asarray(cyt_temp) + np.asarray(stim_vec)    
   moose.vec('/model/chem/dend/cytosol').n = cyt_update




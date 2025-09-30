import moose
import numpy as np
import matplotlib.pyplot as plt
import chemModel as CM
import totalEnergy as TE
import dendShape
from scipy import optimize
from scipy.optimize import LinearConstraint, NonlinearConstraint
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import writeXML
import tag
from multiprocessing import Pool
from functools import partial
import random
from scipy.optimize import Bounds
import magicD as mgd
bounds = Bounds([0.001, 0.01, 0.001], [5.0, 1.5, 5.0]) 


threshold_val = 0.0004
boundary_padding = tag.boundary_padding
disable_gui = True
location_list = []
updated_locations = []
mem_coords = []
num_voxels = len(moose.vec('/model/chem/dend/membrane').n)
linear_constraint1 = LinearConstraint( [[1, 0, 0]], [0.001], [5.0] )
linear_constraint2 = LinearConstraint( [[0, 1, 0]], [0.01], [1.5] )
linear_constraint3 = LinearConstraint( [[0, 0, 1]], [0.001], [5.0] )
n_opt_it = 20

k_agg = tag.k_agg
k_agg2 = tag.k_agg2

def relocate_mem(mem_values, location_relocate):
    profile_indices = np.where(np.asarray(mem_values) > 0)[0]
    profile_values = []
    for pi in range(len(profile_indices)):
        profile_values.append( mem_values[profile_indices[pi]] )
    copy_mem = mem_values.copy()
    for mi in range(len(copy_mem)):
        copy_mem[mi] = 0
    assign_count = 0
    location_relocate_index = int( location_relocate / TE.diffL )
    width = len(profile_indices)
    for ac in  list(np.arange(location_relocate_index - int(0.5 * width), location_relocate_index - int(0.5 * width) + width, 1)):
        copy_mem[ac] = profile_values[assign_count]
        assign_count = assign_count + 1
    return copy_mem
    
    
def branch_select(filename, conc):
   params = pd.read_csv(filename)
   print("Branch parameters: ", params)
   dE = params["p"][0] * conc**2 + params["p"][1]
   print("dE :", dE)
   if dE > 0:
       branch_select_flag = +1
   if dE < 0:
       branch_select_flag = -1
   return branch_select_flag, dE 


for nv in range(num_voxels):
    mem_coords.append(moose.element('/model/chem/dend/mesh['+str(nv)+']').Coordinates[0])


def calc_criterion(conc):
    const1 = (TE.thickness * TE.Na * conc) / TE.phisat
    mu0 = TE.calc_mu0()
    const2 = np.abs(mu0) / TE.RT
    return const1 * np.exp(const2)

def limit_check(ret, phitot, conc, ph_it, Limits, branch):
    if branch == "Low Phi":
        if ret.success and ph_it < Limits[1]:
              print("All conditions satisfied")  
              success_ret = ret
        if not ret.success or ph_it > Limits[1]:
              print("Going to do multi optimization")
              init_samples = calc_failSafe(n_opt_it)
              success_ret = multi_optimization(init_samples, phitot, conc, Limits, branch)
    if branch == "High Phi":
        if ret.success and ph_it > Limits[0]:
              print("All conditions satisfied")  
              success_ret = ret
        if not ret.success or ph_it < Limits[0]:
              print("Going to do multi optimization")
              init_samples = calc_failSafe(n_opt_it)
              success_ret = multi_optimization(init_samples, phitot, conc, Limits, branch)
    return success_ret          


def calc_energy(init_sample, phitot, conc, Limits):    
    print("Init samples: ", init_sample)
    def calc_phi(xp): 
       area = TE.calc_dome_area( (xp[0]/xp[1]) * 1e-6, xp[1])
       phi = phitot / (area * TE.phisat)
       return phi
    nonlinear_constraint = NonlinearConstraint(calc_phi, Limits[0], Limits[1])
    try:
      ret = optimize.minimize( TE.total_energy, init_sample, method='trust-constr', args = (phitot, conc),  constraints = [linear_constraint1, linear_constraint2,linear_constraint3, nonlinear_constraint])
      #ret = optimize.minimize( TE.total_energy, init_sample, method='trust-constr', args = (phitot, conc),  constraints = [linear_constraint1, linear_constraint2,linear_constraint3])
      print("RETURN STATUS: ", ret)
    except:
      print("Contain NaN")
      for i in range(1,n_opt_it):
          try:
             for rs in range(3):
                init_sample[rs] = random.uniform(0.001, 1.0)
             ret = optimize.minimize( TE.total_energy, init_sample, method='trust-constr', args = (phitot, conc),  constraints = [linear_constraint1, linear_constraint2,linear_constraint3, nonlinear_constraint])
             #ret = optimize.minimize( TE.total_energy, init_sample, method='trust-constr', args = (phitot, conc),  constraints = [linear_constraint1, linear_constraint2,linear_constraint3])
             break
          except:
              print("Contain NaN")

    return ret

def calc_failSafe(n_opt):
    init_samples = gen_random(n_opt)
    return init_samples



def gen_random(n_samples):
    total_samples = []
    for n in range(n_samples):
        dx0 = []
        for i in range(3):
              dx0.append(random.uniform(0.001, 1.0))
        total_samples.append(np.asarray(dx0))
    return total_samples

def split_samples(init_samples, num_split):
    samples = []
    start = 0
    each_set_n = int(len(init_samples) / num_split)
    for nsplit in range(num_split):
       sample_set = [] 
       for isamp in range(start, start + each_set_n):
           sample_set.append(init_samples[isamp])
           start = isamp
       samples.append( sample_set )
    return samples   


def multi_optimization(init_samples, phitot, conc, Limits, branch):
    shape_dt, opt_sample_min = calc_convergence()
    print("Doing multi optimization: " + str(opt_sample_min) + " samples at a time") 
    samples = split_samples(init_samples, opt_sample_min)
    for samp in samples:
       print("Split samples: ", samp, len(samp))
       with Pool() as pool:
         result = pool.map(partial(calc_energy, phitot = phitot, conc = conc, Limits = Limits), samp)
       pool.close()
       success_res, success_flag = check_success(result, phitot, Limits, branch)   
       print("Success flag: ", success_flag)
       if success_flag == 1:
           print("Successfull res: ", success_res)
           break
    print("Returned ret from multi opt: ", success_res)
    print("Is it success: ", success_flag)
    return success_res

def check_success(result, phitot, Limits, branch):
    for res in result:
       rm, theta, rp, area, ph_it, alength, r0 = calc_params(res, phitot)
       success_res = res
       success_flag = 0
       if branch == "High Phi":
          if res.success and ph_it > Limits[0]:
             success_res = res
             success_flag = 1
             break
       if branch == "Low Phi":  
          if res.success and ph_it < Limits[1]:
             success_res = res
             success_flag = 1
             break
    return success_res, success_flag

def calc_params(opt_data, phitot):
   print(opt_data) 
   theta = opt_data.x[1]
   alength = opt_data.x[0] * 1e-6
   rm = alength / theta
   rp = opt_data.x[2] * 1e-6
   r0 = (rm + rp) * np.sin(theta)
   area = TE.calc_dome_area(rm, theta)
   ph_it = phitot / (area * TE.phisat)
   return rm, theta, rp, area, ph_it, alength, r0


class essential_lists():
    """
    This class has functions that generate the essential lists.
    I need the lists of starting and ending of domes, and their geometrical parameters

    """
    def reinit(self):
        """
        The constructor
        """
        #Copying the previous phi list
        try:
           self.OldSphi = self.Sphi_list.copy()
        except:
           print("No protrusions")  
           self.OldSphi = []
        try:
           self.oldInitConList = self.initConList.copy()
        except:
           print("No optimization run yet")
           self.oldInitConList = []
        #####################################################
        ## Holds the starting coordinates of the dome (in m)
        self.dome_starts = []
        ## Holds the end coordinates of the dome (in m)
        self.dome_ends = []
        ## Holds the saddle beginnings
        self.saddle_starts = []
        ## Holds the saddle ends (in um)
        self.saddle_ends = []
        ## Holds the starting coordinates of the dome (in voxels)
        self.dome_start_voxels = []
        ## Holds the end coordinates of the dome (in voxels)
        self.dome_end_voxels = []
        ## Holds 1/Rd values of each dome
        self.Hmean_list = []
        self.Rp_list = []
        ## Holds the theta values of each dome
        self.theta_list = []
        ## Holds the sufrace fractions of each dome phi
        self.Sphi_list = []
        ## Holds the area of each dome
        self.area_list = []
        ####################################################
        self.initConList = []

    def __init__(self):
        mem_non_zero_list = list(np.where(moose.vec('/model/chem/dend/membrane').n > 1e-1)[0])
        print("Mem non zero: ", mem_non_zero_list)
        ## List of voxels that are contiguous. Generated from membrane molecule distribution
        self.contiguous_vecs = self.contiguous_chunk(mem_non_zero_list)
        print("CONTIG. VECS")
        print( self.contiguous_vecs )
        self.reinit()

        cytN_zones, cytC_zones = self.content('cytosol', self.contiguous_vecs)
        ## Number of cytosol molecules under each dome
        self.cytN_zones = cytN_zones
        ## Concentration of cytosol molecules under each dome
        self.cytC_zones = cytC_zones
        print("CALCULATED CYT_C: ", self.cytC_zones)
        memN_zones, memC_zones = self.content('membrane', self.contiguous_vecs)
        ## Number of membrane molecules under each dome
        self.memN_zones = memN_zones
        ## Concentration of membrane molecules under each dome
        self.memC_zones = memC_zones
        print("\n ESSENTIAL LISTS CALLED \n")
        num_iterations = 5
        for i in range(num_iterations):
            self.under_dome_vec_list = self.under_dome_vecs()
            print("DOME VECS")
            print( self.under_dome_vec_list )

            self.reinit()

            cytN_zones, cytC_zones = self.content('cytosol', self.under_dome_vec_list)
            ## Number of cytosol molecules under each dome
            self.cytN_zones = cytN_zones
            ## Concentration of cytosol molecules under each dome
            self.cytC_zones = cytC_zones
            print("CALCULATED CYT_C: ", self.cytC_zones)
            memN_zones, memC_zones = self.content('membrane', self.under_dome_vec_list)
            ## Number of membrane molecules under each dome
            self.memN_zones = memN_zones
            ## Concentration of membrane molecules under each dome
            self.memC_zones = memC_zones
            

    def find_saddle_intersection(self, saddle_starts, saddle_ends):
        intersections = [] #spine, intersect pairs
        intersect_coordinates = []
        for sad_l in range(len(saddle_starts) - 1):
            if saddle_ends[sad_l] > saddle_starts[sad_l + 1]:
                print("SADDLE INTERSECTS")
                coordinate_intersection = ( saddle_ends[sad_l] + saddle_starts[sad_l+1] ) / 2.0
                intersect_coordinates.append( coordinate_intersection )
                intersections.append( [sad_l, sad_l + 1] )

        return intersections, intersect_coordinates        

    def merge_spine_geom(self, intersections, intersect_coordinates, N_zones):
        for m in range(len(intersections)):
            update_N = N_zones[ intersections[m][0] ] + N_zones[ intersections[m][1] ]
            update_cytC = (self.cytC_zones[ intersections[m][0] ] + self.cytC_zones[ intersections[m][1] ]) / 2.0
            updated_cyt_conc = update_cytC
            print("MERGED SPINES: ", intersections[m][0], intersections[m][1])
            print("UPDATED MOLECULES: ", update_N)
            #ret = optimize.minimize( TE.total_energy, x0=TE.x0, bounds = (TE.rmbounds,TE.tbounds,TE.rpbounds),args=(update_N, updated_cyt_conc),method = 'L-BFGS-B', options = {'ftol':1e-12} )
            ret = self.random_init_merge(update_N, updated_cyt_conc)

            print("UPDATED GEOMETRY: ", ret.x[0] / ret.x[1], ret.x[1])
            width_new_dome = 2 * (ret.x[0] * 1e-6 / ret.x[1]) * np.sin(ret.x[1])
            dome_start_up = intersect_coordinates[m] - 0.5 * width_new_dome
            dome_end_up = intersect_coordinates[m] + 0.5 * width_new_dome
            Rs = ret.x[2] * 1e-6
            width_new_saddle = 2 * Rs * np.sin(ret.x[1])
            saddle_start_up = intersect_coordinates[m] - 0.5 * width_new_dome - 0.5 * width_new_saddle
            saddle_end_up = intersect_coordinates[m] + 0.5 * width_new_dome + 0.5 * width_new_saddle

            print("UPDATED DOME COORDS: ", dome_start_up, dome_end_up)    
            self.saddle_starts[ intersections[m][0] ] = saddle_start_up
            self.saddle_ends[ intersections[m][0] ] = saddle_end_up
            self.dome_starts[ intersections[m][0] ] = dome_start_up
            self.dome_ends[ intersections[m][0] ] = dome_end_up
            self.dome_start_voxels[ intersections[m][0] ] = int( dome_start_up / TE.diffL )
            self.dome_end_voxels[ intersections[m][0] ] = int( dome_end_up / TE.diffL )
            location_list[ intersections[m][0] ] = ( dome_start_up + dome_end_up ) / 2.0
            self.dome_starts.remove( self.dome_starts[ intersections[m][1] ] )
            self.dome_ends.remove( self.dome_ends[ intersections[m][1] ] )
            self.dome_start_voxels.remove( self.dome_start_voxels[ intersections[m][1] ] )
            self.dome_end_voxels.remove( self.dome_end_voxels[ intersections[m][1] ] )
            self.saddle_starts.remove( self.saddle_starts[ intersections[m][1] ] )
            self.saddle_ends.remove( self.saddle_ends[ intersections[m][1] ] )
            try:
                  location_list.remove( location_list[ intersections[m][1] ] )
            except:
                  print("No index for merge found") 
            print("UPDATED COORDINATES")
            print("INTERSECTION: ", intersections)
            print("SADDLE STARTS: ", self.saddle_starts)
            print("DOME STARTS: ", self.dome_starts)
            print("DOME ENDS: ", self.dome_ends)
            print("SADDLE ENDS: ", self.saddle_ends)
    
    def random_init_merge(self, updated_phitot, updated_cytConc):
       init_guess = TE.x0
       phitot = updated_phitot
       conc = updated_cytConc
       branch_select_flag = -1

       if branch_select_flag > 0:
         phi_Llim = 0.0
         phi_Hlim = 0.1
         branch = "Low Phi"
       if branch_select_flag < 0:
         phi_Llim = 0.2
         phi_Hlim = 1.0
         branch = "High Phi"
       Limits = [phi_Llim, phi_Hlim]
       print("Limits: ", Limits)
       print("Entering optimization")
       ret = calc_energy(init_guess, phitot, conc, Limits)
       rm, theta, rp, area, ph_it, alength, r0 = calc_params(ret, phitot)
       success_ret = limit_check(ret, phitot, conc, ph_it, Limits, branch)
       return success_ret

    def random_init(self, phitot_plus_indices):
       init_guess = TE.x0
       cyt_concs = self.cytC_zones
       print(phitot_plus_indices)
       phitot = phitot_plus_indices[0]
       conc = cyt_concs[phitot_plus_indices[1]]
       branch_select_flag = -1

       if branch_select_flag > 0:
         phi_Llim = 0.0
         phi_Hlim = 0.1
         branch = "Low Phi"
       if branch_select_flag < 0:
         phi_Llim = 0.2
         phi_Hlim = 1.0
         branch = "High Phi"
       Limits = [phi_Llim, phi_Hlim]
       print("Limits: ", Limits)
       print("Entering optimization")
       ret = calc_energy(init_guess, phitot, conc, Limits)
       rm, theta, rp, area, ph_it, alength, r0 = calc_params(ret, phitot)
       success_ret = limit_check(ret, phitot, conc, ph_it, Limits, branch)
       return success_ret

  

    def content(self, pool_name, summation_vecs):
        """
        This is the function that counts and assigns the number and concentrations of molecules.
        pool_name is either cytosol or membrane.
        """
        N_zones = []
        C_zones = []
        location_list = []
        mem_count = []
        temp_pool = moose.vec('/model/chem/dend/' + pool_name).n
        temp_pool_conc = moose.vec('/model/chem/dend/' + pool_name).conc
        for i in range(len(summation_vecs)):
            temp_count = sum(temp_pool[summation_vecs[i]])
            length_of_zones = len(summation_vecs[i])
            if pool_name == 'membrane':
                 total_local_conc = (self.cytN_zones[i] + temp_count) / (np.pi * TE.Rc**2 * TE.diffL * length_of_zones * TE.Na)
                 print("LOCAL CONC" + str(i) +": ", 1e3 * total_local_conc)
                 if current_time > tag.spine_create_delay:
                    if total_local_conc > threshold_val:
                       mem_count.append([temp_count, i])
                 else:
                    mem_count.append([temp_count, i])
            temp_count_conc = sum(temp_pool_conc[summation_vecs[i]])
            N_zones.append(temp_count)
            print(pool_name + " NZONES: ", N_zones)
            if pool_name == 'cytosol':
                C_zones.append( temp_count_conc / length_of_zones )
        if pool_name == 'membrane':
               rets = []
               print("Length of mem Zones: ", len(mem_count))
               with ThreadPoolExecutor(max_workers = tag.num_spines + 2) as executor:
                   ret_values = executor.map(self.random_init, mem_count)

               for rt in ret_values:
                   rets.append(rt)
               print("Length of rets: ", len(rets))
               for rt in range(len(mem_count)):  
                   print("Spine index: ", mem_count[rt][1])
                   area = 2 * np.pi * (rets[rt].x[0] * 1e-6 / rets[rt].x[1])**2 * (1 - np.cos(rets[rt].x[1]))
                   self.area_list.append(area)
                   self.Sphi_list.append( mem_count[rt][0] / (area * TE.phisat) )
                   self.Hmean_list.append( 1 / (rets[rt].x[0] * 1e-6 / rets[rt].x[1]) )
                   self.Rp_list.append( rets[rt].x[2] * 1e-6 )
                   self.theta_list.append( rets[rt].x[1] )
                   volume = area * TE.thickness
                   C_zones.append( mem_count[rt][0] / (TE.Na * volume) )
                   width_new_dome = 2 * (rets[rt].x[0] * 1e-6 / rets[rt].x[1]) * np.sin(rets[rt].x[1])
                   Rs = rets[rt].x[2] * 1e-6
                   width_new_saddle = 2 * Rs * np.sin(rets[rt].x[1])
                # Finding the center of protrusion
                   center_voxel = int( len(summation_vecs[mem_count[rt][1]]) / 2.0)
                   center_voxel_in_um = summation_vecs[mem_count[rt][1]][center_voxel] * TE.diffL
                # Storing the location in location_list
                   location_list.append( center_voxel_in_um )
                   self.dome_starts.append( center_voxel_in_um - 0.5 * width_new_dome) 
                   self.saddle_starts.append( center_voxel_in_um - 0.5 * width_new_dome - 0.5 * width_new_saddle )
                   self.dome_ends.append( center_voxel_in_um + 0.5 * width_new_dome ) 
                   self.saddle_ends.append( center_voxel_in_um + 0.5 * width_new_dome + 0.5 * width_new_saddle )
                   self.dome_start_voxels.append( int( (center_voxel_in_um - 0.5 * width_new_dome ) / TE.diffL ) )
                   self.dome_end_voxels.append( int( (center_voxel_in_um + 0.5 * width_new_dome ) / TE.diffL ) )
        if pool_name == 'membrane':
               intersections, intersect_coordinates = self.find_saddle_intersection(self.saddle_starts, self.saddle_ends)        
               print("INTERSECTION: ", intersections)
               self.merge_spine_geom(intersections, intersect_coordinates, N_zones)
        return N_zones, C_zones
    
    def contiguous_chunk(self,list_item):
        """
        This is the function that looks for contiguous chunks of voxels
        that has non zero membrane molecules.
        list_item has the nonzero membrane voxels of the dendrite
        """
        contiguous_vecs = []
        temp_vec = []
        augmented_list = list_item + [0]
        for i in range(len(list_item)):
            if (augmented_list[i + 1] == augmented_list[i] + 1):
                temp_vec.append(augmented_list[i])
            else:
                temp_vec.append(augmented_list[i])
                contiguous_vecs.append(temp_vec)
                temp_vec = []
        return contiguous_vecs    
    
    def under_dome_vecs(self):
        under_dome_vecs = []
        for i in range(len(self.dome_start_voxels)):
            voxels_under_dome = list( np.arange(self.dome_start_voxels[i], self.dome_end_voxels[i]) )
            voxels_under_dome.append( self.dome_end_voxels[i] )
            under_dome_vecs.append(voxels_under_dome)
        return under_dome_vecs     
            





class update():
    """
    Calling this class updates
    1) Geometrical field - Hmean, theta
    2) Concentration field for calculation of Nernst potential - memConc and cytConc
    3) Surface fraction phi - Sphi
    4) Updates the marker of the dome
    All these fields are assigned to the MOOSE pools that holds them
    """
    def __init__(self,dome_start_v, dome_end_v, Hmean_v, theta_v, Sphi_v, memConc_v, cytConc_v, memN_zones_v):
        """
        Constructor that does the updating. Uses the lists created by essential_lists
        """
        marker = np.zeros(len(moose.vec('/model/chem/dend/marker').n))
        Hmean = np.zeros(len(moose.vec('/model/chem/dend/Hmean').n))
        theta = np.zeros(len(moose.vec('/model/chem/dend/theta').n))
        Sphi = np.zeros(len(moose.vec('/model/chem/dend/Sphi').n))
        memConc = np.zeros(len(moose.vec('/model/chem/dend/memConc').n))
        cytConc = np.zeros(len(moose.vec('/model/chem/dend/cytConc').n))
        memv = np.zeros(len(moose.vec('/model/chem/dend/memv').n))
        for d in range(len(dome_start_v)):
            print("START VOXEL: ", dome_start_v[d])
            print("END VOXELS: ", dome_end_v[d])
            marker[dome_start_v[d] : dome_end_v[d] + 1] = 1.0
            Hmean[dome_start_v[d] : dome_end_v[d] + 1] = Hmean_v[d]
            theta[dome_start_v[d] : dome_end_v[d] + 1] = theta_v[d]
            Sphi[dome_start_v[d] : dome_end_v[d] + 1] = Sphi_v[d]
            memConc[dome_start_v[d] : dome_end_v[d] + 1] = memConc_v[d]
            cytConc[dome_start_v[d] : dome_end_v[d] + 1] = cytConc_v[d]
            memv[dome_start_v[d] : dome_end_v[d] + 1] = memN_zones_v[d]
        moose.vec('/model/chem/dend/marker').n = marker
        moose.vec('/model/chem/dend/Hmean').n = Hmean
        moose.vec('/model/chem/dend/theta').n = theta
        moose.vec('/model/chem/dend/Sphi').n = Sphi
        moose.vec('/model/chem/dend/memConc').n = memConc
        moose.vec('/model/chem/dend/cytConc').n = cytConc
        moose.vec('/model/chem/dend/memv').n = memv

                
    


class new_spine_attrib(): 
    """
    This class has functions and data to generate a new spine. 
    If the spine is an inactive dummy, then width_um parameter has to be supplied. 
    phitot = 1 to prevent null phi. 
    If not a dummy, then phitot has to be the desired value,
    i.e., the number of molecules in the dome.
    """
    #When a spine is initiated there is no molecules in it. So I am not initializing memConc here. To prevent the chemical potential from acting, I am not initializing cytConc either. log(1e-24 / 1e-24) = log(1) = 0
    def __init__(self, location_um, width_um, theta, dummy):
         """
         This is the constructor. Generates a dummy or non dummy spine parameters,
         i.e., theta, Rd, phi. Marks the region of dome.
         Spine coordinates are calculated using dendShape from dendShape.py.
         Distributes Hmean, theta, Sphi fields by calling init_geom_par()
         dummy[True] = 0 means not a dummy
         dummy[True] = 1 means a dummy
         If dummy, membrane molecules are distributed by init_mem()
         """
         print("Length in voxels: ", int(width_um / TE.diffL))
         location_list.append(location_um)
         location_list.sort()
         if dummy[True] == 1:
            ## the width of a dummy spine (in m) 
            self.width_um = width_um 
            ##The number of molecules inside the dome (in numbers)
            self.phitot = dummy['phitot']
            ##The angle that tangent to the inflexion point makes with the horizontal.
            self.theta = theta
            rm_p_rp = self.width_um / np.sin(theta)
            self.rm = (rm_p_rp / 2.0) * 1e6
            self.rp = (rm_p_rp / 2.0) * 1e6
            ## The surface fraction = phitot / (area_dome * phisat)
            self.phi = self.phitot / ( TE.phisat * 2 * np.pi * (self.rm * 1e-6)**2 * (1 - np.cos(self.theta)))

         self.x_vec = np.linspace(0, TE.Length, num_voxels)
         ## The starting coordinate of the spine (in voxel number)
         self.start_voxel = int( ( location_um -  self.width_um / 2.0) / CM.diffL )
         ## The ending coordinate of the spine (in voxel number)
         self.end_voxel = int( ( location_um +  self.width_um / 2.0) / CM.diffL )
         ## Just a temperory variable used to assign value = 1 for voxels under the dome
         self.Mvec = np.asarray( [0] * num_voxels )
         self.Mvec[self.start_voxel : self.end_voxel + 1] = 1
         moose.vec('/model/chem/dend/marker').n = moose.vec('/model/chem/dend/marker').n + self.Mvec
         x, y, dome_start, dome_end = dendShape.dendShape( location_um, self.theta, self. rm, self.rp )
         ## dome_start is the start coordinate of dome (in um)
         self.dome_start = dome_start
         ## dome_end is the end coordinate of dome (in um)
         self.dome_end = dome_end
         print("initial corrds: ", dome_start, dome_end)
         ## x has the x-coordinates of the spine (in m)
         self.x = x[np.where(x <= 10)]
         ## y has the y-coordinates of the spine (in m)
         self.y = y[np.where(x <= 10)] 
         self.init_geom_par()

         if dummy[True] == 0:
             self.init_mem()



    def init_geom_par(self):
         """
         The chemical model need to know certain parameters for activating Kf,
         i.e., recruitment rate from cytosol to membrane. 
         These parameters should be uniformly distributed in the dome region.
         That is what this functions does
         """
         thetaVec = moose.vec( '/model/chem/dend/theta' ).n
         HmeanVec = moose.vec( '/model/chem/dend/Hmean' ).n
         SphiVec = moose.vec( '/model/chem/dend/Sphi' ).n
         for i in range(num_voxels):
           thetaVec[i] = self.Mvec[i] * self.theta
           HmeanVec[i] = self.Mvec[i] * ( 1 / (self.rm * 1e-6) )
           SphiVec[i] = self.Mvec[i] * self.phi
         moose.vec( '/model/chem/dend/theta' ).n = moose.vec( '/model/chem/dend/theta' ).n + thetaVec
         moose.vec( '/model/chem/dend/Hmean' ).n = moose.vec( '/model/chem/dend/Hmean' ).n + HmeanVec
         moose.vec( '/model/chem/dend/Sphi' ).n = moose.vec( '/model/chem/dend/Sphi' ).n + SphiVec

    def init_mem(self):
        """
        If the spine is not a dummy, membrane molecules are present in it.
        Then we need to distribute them accordingly, i.e., uniformly inside the dome.
        Each voxel has phitot / length_of_protrusion
        """
        membrane_old = list(moose.vec( '/model/chem/dend/membrane' ).n)
        update_mem = membrane_old.copy()
        for i in range(num_voxels):
            update_mem[i] = self.Mvec[i] * self.phitot / (self.end_voxel - self.start_voxel + 1)
        moose.vec( '/model/chem/dend/membrane' ).n = np.asarray(update_mem) + np.asarray(membrane_old)    
        print("Membrane sum after new spine: ", sum(moose.vec( '/model/chem/dend/membrane' ).n))
        print("CYtosol sum after new spine: ", sum(moose.vec( '/model/chem/dend/cytosol' ).n))


    
    def plot_mem(self):
         plt.plot(self.x_vec, moose.vec('/model/chem/dend/marker').n)
         plt.show()

    def plot_shape(self):
         plt.plot(self.x, self.y)
         plt.show()


fig, ( (ax1, ax2) , (ax3, ax4), (ax5, ax6), (ax7, ax8) ) = plt.subplots(4 , 2)

ax1.set_title("Hmean")
ax2.set_title("Membrane")
ax3.set_title("Cytosol")
ax4.set_title("$\phi$")
ax5.set_title("memConc")
ax6.set_title("cytConc")
ax7.set_title("Shapes")
ax8.set_title("Dome number")

runtime = tag.runtime
phi_save = []
shapeTimes = []
Hmean_save = []
theta_save = []
mConc_spine_wise = []
cConc_spine_wise = []
location_save = []
cyt_save = []
cyt_conc_save = []
mem_save = []
marker_save = []
hmean_save = []
Sphi_save = []
mConc_save = []
cConc_save = []
kf_save = []
pool_kf_save = []
pot_save = []
saddle_start_save = []
Rp_save = []
Rac_save = []
Tiam_save = []
Ca_save = []
dimer_save = []
thr286_save = []



num_spines = tag.num_spines

def calc_effective_conc(cytC, memN, dome_start_voxel, dome_end_voxel):
    total_conc_dome = []
    total_conc_dend = []
    for i in range(len(cytC)):
       num_voxels = int(dome_end_voxel[i] - dome_start_voxel[i])
       volume = np.pi * (0.5 * TE.dendDia)**2 * num_voxels * TE.diffL
       cytN = cytC[i] * volume * TE.Na
       total_molecules_dome = cytN + memN[i]
       total_conc_dome.append(1e3 * total_molecules_dome / (TE.Na * volume)) #in um
       sum_cytosol = sum(moose.vec('/model/chem/dend/cytosol').n)
       total_conc_dend.append(1e3 * (memN[i] + sum_cytosol) / (TE.Na * TE.V_cyt))
    return total_conc_dome, total_conc_dend

def calc_volume_voxel(num_voxels):
    voxels_volume = []
    for nv in range(num_voxels):
         voxels_volume.append( moose.vec('/model/chem/dend/mesh')[nv].volume )
    print("Voxel volumes: ", voxels_volume)     
    return voxels_volume 


def threshold_detect(loc, threshold_val):
    cytosol = moose.vec('/model/chem/dend/cytosol').conc
    loc_voxel = int( loc / TE.diffL )
    if cytosol[loc_voxel] > threshold_val:
        flag = 1
    else:
        flag = 0
    return flag, cytosol[loc_voxel]     


print("Dimer initial: ", moose.vec('/model/chem/dend/IRSp53_dimer').n)


def remove_duplicates(locations):
    res = []
    for i in locations:
        if i not in res:
            res.append(i)
    return res        


if tag.chemTurnOff == False:
   glut_stim_times = tag.glut_stim_times
   glut_stim_locations = tag.glut_stim_locations
   df_stim = pd.DataFrame()
   df_stim["stimTimes"] = glut_stim_times
   df_stim["stimLocations"] = glut_stim_locations
   df_stim.to_csv("./" + tag.tag_string + "stimD.csv")
inter_spine_distance = 0.7e-6
dummy = {True : 1, 'phitot' : 1.0}  #Careful to change the phitot when dummy[True] = False because then it is a mature spine with molecules inside it.



def check_for_max():
      conc_list = list( moose.vec('/model/chem/dend/cytosol').conc )
      left_boundary_voxel = int(boundary_padding / TE.diffL)
      right_boundary_voxel = int( (tag.Length * 1e-6 - boundary_padding) / TE.diffL )
      conc_list = conc_list[left_boundary_voxel : right_boundary_voxel]
      index_max = left_boundary_voxel + conc_list.index( max(conc_list) )
      max_conc = max( moose.vec('/model/chem/dend/cytosol').conc )
      print("Reached max loc: ", index_max * TE.diffL)
      print("Max conc: ", max_conc)
      return index_max * TE.diffL, max_conc       

def check_in_domain(location):
    search_limit = boundary_padding
    if location > search_limit and location < tag.Length * 1e-6 - search_limit:
        dflag = 1
    else:
        dflag = 0
    return dflag    

def check_proximity(location, tolerance):
    pr_flag = -1
    print("Checking proximity: ", location)
    print("UPDATED LOCATIONS in check_proximity: ", updated_locations)
    for u_l in updated_locations:
        if abs(location - u_l) < tolerance:
           pr_flag = 0
           break
    return pr_flag

dummy_locations = tag.dummy_spine_locations
spine_create_times = tag.spine_create_times

def init_dummy_create(n):
      print("Creating spine")
      new_spine_attrib(dummy_locations[n], 0.5e-6, 0.2, dummy)
      updated_locations.append(dummy_locations[n])

def stimless_dummy(n):
  if tag.init_dummy:
     init_dummy_create(0)
     dummy_locations.pop(0)
     spine_create_times.pop(0)


dt_list = []
if tag.chemTurnOff == False:
   stim_start_time = 0
   stim_end_time = tag.bg_end   
else:
   stim_start_time = 0
   stim_end_time = tag.runtime
current_time = 0
#for r in range( int(runtime/shape_dt) ):

if tag.chemTurnOff == False:
  deletion_delay = tag.deletion_delay
  try:
    deletion_time = glut_stim_times[-1] + deletion_delay
  except:
    deletion_time = 0
else:
  deletion_delay = 0
  deletion_time = 0
rand_prob_list = []
rand_prob_spines = []

def calc_convergence():
   voxel_volume = np.pi * (0.5 * TE.dendDia)**2 * TE.diffL
   max_cytosol = 1e3 * max(moose.vec('/model/chem/dend/cytosol').n) / (TE.Na * voxel_volume)
   min_cytosol = 1e3 * min(moose.vec('/model/chem/dend/cytosol').n) / (TE.Na * voxel_volume)
   diff_conc_cytosol = max_cytosol - min_cytosol
   print("Min cyt: ", min_cytosol)
   print("Max cyt: ", max_cytosol)
   print("diff cyt: ", diff_conc_cytosol)
   if diff_conc_cytosol < 1e3 * tag.cyt_conc:
         shape_dt = 0.2
         opt_sample_min = 2
   else:
         shape_dt = 0.6
         opt_sample_min = 10
   return shape_dt, opt_sample_min
      

deleted_spine = []
new_spine_created_time = 1.0 
while current_time < runtime:
   opt_sample_min = 10
   if current_time > stim_start_time and current_time < stim_end_time:
      #shape_dt = tag.CaTau
      shape_dt = 0.2
      if tag.chemTurnOff == True:
         shape_dt, opt_sample_min = calc_convergence()
   else:
      shape_dt, opt_sample_min = calc_convergence()
   print("Shape dt: ", shape_dt)
   dt_list.append(current_time)
   sum_cytosol = sum(moose.vec('/model/chem/dend/cytosol').n)
   cytosol_conc = sum_cytosol / (TE.V_cyt * TE.Na)
   if len(dummy_locations) > 0:
      if cytosol_conc > threshold_val:
         for ns in range(tag.num_spines):
            stimless_dummy(0)
      
   cytosol_conc_array = moose.vec('/model/chem/dend/cytosol').conc

   def spec_loc_create(old_spine_created_time):
       spine_created_time = old_spine_created_time
       if len(dummy_locations) > 0:
          if cytosol_conc > threshold_val and current_time > spine_create_times[0]:
             print("SPINE CREATION CONDITIONS EXIST")
             print("CREATING SPINE AT: ", dummy_locations[0])
             init_dummy_create(0)
             dummy_locations.pop(0)
             spine_create_times.pop(0)
             spine_created_time = current_time
       return spine_created_time     
  

   def search_opt_loc(old_spine_created_time):
     spine_created_time = old_spine_created_time
     max_loc, max_conc = check_for_max()
     print("Max location and value: ", max_loc, max_conc)
     if len(updated_locations) > 0:
          pr_flag = check_proximity(max_loc,  0.5e-6)
     else:
          pr_flag = 1
     print("Proximity flag: ", pr_flag)
     dflag = check_in_domain(max_loc)
     print("In domain flag: ", dflag)
     if dflag == 1:
        if pr_flag == 1 or pr_flag == -1:
           print("CREATING SPINE")
           new_spine_attrib(max_loc, 0.5e-6, 0.2, dummy)
           spine_created_time = current_time
           updated_locations.append(max_loc)
           #break
        else:
           print("SPINE EXIST AT THIS LOCATION: ", updated_locations)
     else:          
       print("SPINE OUTSIDE DOMAIN: ", updated_locations)
       print("Exiting the spine creatng loop")
     return spine_created_time

   if tag.chemTurnOff == False:
      print("Old spine creation time: ", new_spine_created_time)
      if (current_time - new_spine_created_time) >= shape_dt:
          print("Searching for threshold: current time: previous time: ", current_time, new_spine_created_time)
          print("CYTOSOL AND THRESHOLD: ", cytosol_conc, threshold_val)
          if cytosol_conc > threshold_val:
             #new_spine_created_time = spec_loc_create(new_spine_created_time)            
             new_spine_created_time = search_opt_loc(new_spine_created_time)
 
             print("New spine creation time: ", new_spine_created_time)
   print(moose.vec('/model/chem/dend/membrane').n)

   print("MEMBRANE SUM: ", sum(moose.vec('/model/chem/dend/membrane').n))
   print("CYTOSOL SUM: ", sum(moose.vec('/model/chem/dend/cytosol').n))
   print("INACTIVE SUM: ", sum(moose.vec('/model/chem/dend/IRSp53_dimer').n))
   print("Ca Max: ", np.max(moose.vec('/model/chem/dend/Ca').n))
   moose.start(shape_dt)


   print("CURRENT TIME: ", current_time)
   print("LOCATION LIST: ", location_list)
   print("UPDATED LOCATION LIST: ", updated_locations)
 
   if len(updated_locations) > 0:
      #if current_time % shape_dt == 0: 
         s = essential_lists()
         mgd.magic_diff(s.dome_start_voxels, s.dome_end_voxels)
         if current_time >= deletion_time and tag.allowDeletion:
               effective_concs_dome, effective_concs_dend = calc_effective_conc(s.cytC_zones, s.memN_zones, s.dome_start_voxels, s.dome_end_voxels)
               print("Effective concs: ", effective_concs_dend)
               location_under_threshold =  np.where(np.asarray(effective_concs_dend) < threshold_val * 1e3)[0]
               if len(location_under_threshold) > 0:
                  print("Locations under threshold: ", location_under_threshold)
                  if tag.tiam_dependency:
                    concs_Tiam_frac = []
                    concs_Tiam = []
                    Tiam_vec = moose.vec('/model/chem/dend/Tiam1').n
                    for scrd in range(len(s.dome_start_voxels)):
                      voxel_nums = list( np.arange(s.dome_start_voxels[scrd], s.dome_end_voxels[scrd], 1) )
                      num_mols_Tiam = sum(Tiam_vec[voxel_nums])
                      num_voxels_Tiam = len(voxel_nums)
                      concs_Tiam.append(1e3 * num_mols_Tiam / (TE.Na * np.pi * TE.Rc**2 * TE.diffL * num_voxels_Tiam))
                      print("Conc Tiam: ", concs_Tiam)
                    concs_Tiam_frac = np.asarray(concs_Tiam) / sum(concs_Tiam)
                    spine_for_deletion = np.random.choice(location_under_threshold)
                    print("Spine for deletion: ", spine_for_deletion)
                    random_elimination_probability = np.random.random()
                    print("Random probability: ", random_elimination_probability)
                    print("Underlying Tiam conc: ", concs_Tiam[spine_for_deletion])
                    if random_elimination_probability > concs_Tiam_frac[spine_for_deletion]: 
                       delete_voxels = list(np.arange(s.dome_start_voxels[spine_for_deletion] - 10, s.dome_end_voxels[spine_for_deletion] + 10,1))
                       print("Deleting spine: ", spine_for_deletion)
                       deleted_spine.append(updated_locations[spine_for_deletion])
                       mgd.delete_spine(delete_voxels)
                       s = essential_lists()
                       u1 = np.random.random()
                       deletion_time = current_time + (1 / sum(concs_Tiam_frac)) * np.log(1 / u1)
                       print("Next deletion time: ", deletion_time)
                           
         updated_locations = []
         for i in range(len(s.dome_starts)):
             updated_locations.append( ( s.dome_starts[i] + s.dome_ends[i] ) / 2.0 )
         ks = moose.element( '/model/chem/dend/ksolve' )
         #rates2 = np.asarray(ks.rateVec['/model/chem/dend/reacF'])
         shapeTimes.append(current_time)
         phi_save.append(s.Sphi_list)
         Hmean_save.append(s.Hmean_list)
         Rp_save.append(s.Rp_list)
         saddle_start_save.append(s.saddle_starts)
         theta_save.append(s.theta_list)
         mConc_spine_wise.append(s.memC_zones)
         cConc_spine_wise.append(s.cytC_zones)
         print("LOCATION SAVE: ", location_save)
         print("LOCATION LIST: ", location_list) 
         location_save.append(updated_locations)
         update(s.dome_start_voxels, s.dome_end_voxels, s.Hmean_list, s.theta_list, s.Sphi_list, s.memC_zones, s.cytC_zones, s.memN_zones)
         print("Molecules in dome: ", s.memN_zones)
         print("Molecules under dome: ", s.cytN_zones)
         print("Phi list: ", s.Sphi_list)
         dphi = pd.DataFrame(phi_save)
         dHmean = pd.DataFrame(Hmean_save)
         dRp = pd.DataFrame(Rp_save)
         dtheta = pd.DataFrame(theta_save)
         dmConc = pd.DataFrame(mConc_spine_wise)
         dcConc = pd.DataFrame(cConc_spine_wise)
         dlocations = pd.DataFrame(location_save)
         dss = pd.DataFrame(saddle_start_save)
            
   Rac_save.append(moose.vec('/model/chem/dend/Rac_GTP').n)
   Tiam_save.append(moose.vec('/model/chem/dend/Tiam1').n)
   thr286_save.append(moose.vec('/model/chem/dend/CaMKII_gr/tot_CaM_CaMKII').n)
   Ca_save.append(moose.vec('/model/chem/dend/Ca').n)
   cyt_save.append(moose.vec('/model/chem/dend/cytosol').n)
   dimer_save.append(moose.vec('/model/chem/dend/IRSp53_dimer').n)
   cyt_conc_save.append(moose.vec('/model/chem/dend/cytosol').conc)
   mem_save.append(moose.vec('/model/chem/dend/membrane').n)
   marker_save.append(moose.vec('/model/chem/dend/marker').n)
   hmean_save.append(moose.vec('/model/chem/dend/Hmean').n)
   Sphi_save.append(moose.vec('/model/chem/dend/Sphi').n)
   mConc_save.append(moose.vec('/model/chem/dend/memConc').n)
   cConc_save.append(moose.vec('/model/chem/dend/cytConc').n)
   current_time = current_time + shape_dt
             #pot_save.append(moose.vec('/model/chem/dend/pot_eval').n)
             #kf_save.append(rates2)
             #pool_kf_save.append(moose.vec('/model/chem/dend/saved_kf').n)

tag_string = tag.tag_string
df_time = pd.DataFrame()
df_time["time"] = dt_list
df_time.to_csv(tag_string + "chemTimes.csv")
df_rand = pd.DataFrame()
df_rand["rand"] = rand_prob_list
df_rand["spine"] = rand_prob_spines
df_rand.to_csv("./" + tag_string + "rand.csv")

ddeleted = pd.DataFrame(deleted_spine)
dphi = pd.DataFrame(phi_save)
dHmean = pd.DataFrame(Hmean_save)
dRp = pd.DataFrame(Rp_save)
dtheta = pd.DataFrame(theta_save)
dmConc = pd.DataFrame(mConc_spine_wise)
dcConc = pd.DataFrame(cConc_spine_wise)
dlocations = pd.DataFrame(location_save)
dss = pd.DataFrame(saddle_start_save)
dtime = pd.DataFrame(shapeTimes)

try:
  dphi.to_csv(tag_string + "phi.csv")
  ddeleted.to_csv(tag_string + "deleted.csv")
  dtime.to_csv(tag_string + "shapeTimes.csv")

  dRp.to_csv(tag_string + "rp.csv")
  #dHmean.columns = labels
  dHmean.to_csv(tag_string + "hmean.csv")

  dss.to_csv(tag_string + "sstart.csv")

  #dtheta.columns = labels
  dtheta.to_csv(tag_string + "theta.csv")

  #dmConc.columns = labels
  dmConc.to_csv(tag_string + "mConc.csv")
      
  #dcConc.columns = labels
  dcConc.to_csv(tag_string + "cConc.csv")

  dlocations.to_csv(tag_string + "locations.csv")
except:
    print("spine not formed")
     
writeXML.writeXML_td(cyt_save, tag_string + "cytosol.xml")
writeXML.writeXML_td(cyt_conc_save, tag_string + "cytosolConc.xml")
writeXML.writeXML_td(mem_save, tag_string + "membrane.xml")
writeXML.writeXML_td(marker_save, tag_string + "marker.xml")
writeXML.writeXML_td(Sphi_save, tag_string + "phi.xml")
writeXML.writeXML_td(hmean_save, tag_string + "hmean.xml")
writeXML.writeXML_td(mConc_save, tag_string + "mConc.xml")
writeXML.writeXML_td(cConc_save, tag_string + "cConc.xml")
writeXML.writeXML_td(kf_save, tag_string + "kf.xml")
writeXML.writeXML_td(pool_kf_save, tag_string + "pool_kf.xml")
writeXML.writeXML_td(pot_save, tag_string + "pot.xml")
writeXML.writeXML_td(Rac_save, tag_string + "Rac.xml")
writeXML.writeXML_td(dimer_save, tag_string + "dimer.xml")
writeXML.writeXML_td(Ca_save, tag_string + "Ca.xml")
writeXML.writeXML_td(thr286_save, tag_string + "thr286.xml")
writeXML.writeXML_td(Tiam_save, tag_string + "Tiam.xml")




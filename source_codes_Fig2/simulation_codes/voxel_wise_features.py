import moose
import pandas as pd
from scipy import optimize
import totalEnergy as TE
import numpy as np

class voxel_wise_features():
    """
    This class has functions that returns the voxel-wise properties
    1) geometrical properties  -  Rd, theta
    2) Concentrations - cytosol and membrane
    3) phi
    """
    def __init__(self, contiguous_vecs, dome_start_voxels, dome_end_voxels, current_time):
        """
        Here, I initialize the afore-mentioned properties.

        Each element in membrane_dist_list is another list with
        the distribution of the number of molecules under the dome.

         
        """
        ## Takes time from MOOSE
        self.current_time = current_time
        ## Taken from makeNewSpine. dome_start position in terms of voxels
        self.dome_start_voxels = dome_start_voxels
        ## Taken from makeNewSpine. dome_end position in terms of voxels
        self.dome_end_voxels = dome_end_voxels
        ## To store areas of each voxels under all domes. It is a list of lists. Each list corresponds to one of the domes
        self.area_voxel_all = []
        ## To store areas of each voxels under the dome.
        self.area_voxel = []
        ## Stores the values of all voxels with non-zero membrane molecules
        self.contiguous_vecs = contiguous_vecs
        ## List of lists. Each list corresponds to the Hmean distribution under a dome
        self.Rd_voxel_all = []
        self.Rd_voxel = []
        ## List of lists. Each list corresponds to the theta distribution under a dome
        self.theta_voxel_all = []
        self.theta_voxel = []
        ## List of lists. Each list corresponds to the phi distribution under a dome.
        self.Sphi_voxel_all = []
        self.Sphi_voxel = []

        self.voxelwise_content()

    def voxelwise_content(self):
        """
        contiguous_vecs is enough for this. Optimization to be run for
        voxelwise content.
        """

        mem_pool = moose.vec('/model/chem/dend/membrane').n
        for cv in range(len(self.contiguous_vecs)):
            self.Rd_voxel = []
            time_list = []
            for cvi in range(len(self.contiguous_vecs[cv])):
               if self.contiguous_vecs[cv][cvi] > self.dome_start_voxels[cv] and self.contiguous_vecs[cv][cvi] < self.dome_end_voxels[cv]:
                 mem_voxel_content = mem_pool[self.contiguous_vecs[cv][cvi]]
                 ret = optimize.minimize( TE.total_energy, x0=TE.x0, bounds = (TE.rmbounds,TE.tbounds),args=(mem_voxel_content,TE.phisat, TE.khat, TE.fs, TE.k, 1, 0, 0, 0, 0.1e-3),method = 'L-BFGS-B', options = {'ftol':1e-12} )
                 self.Rd_voxel.append( ret.x[0] )
                 self.theta_voxel.append( ret.x[1] )
                 area_voxel = 2 * np.pi * ( ret.x[0] * 1e-6 )**2 * ( 1 - np.cos(ret.x[1]) ) 
                 self.area_voxel.append( area_voxel )
                 phi_voxel = mem_voxel_content / ( area_voxel * TE.phisat )
                 self.Sphi_voxel.append( phi_voxel )
                 time_list.append(self.current_time)
            df = pd.DataFrame()
            df["Rd"] = self.Rd_voxel
            df["time"] = time_list
            df.to_csv('./Rd' + str(cv) + str(self.current_time) + '.csv', mode = 'a')


            self.area_voxel_all.append(self.area_voxel)   
            self.Rd_voxel_all.append(self.Rd_voxel)
            self.theta_voxel_all.append(self.theta_voxel)


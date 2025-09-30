import numpy as np
import moose


def magic_diff(dome_start_voxels, dome_end_voxels):
        for cv in range(len(dome_start_voxels)):
           membrane_old = moose.vec('/model/chem/dend/membrane').n
           num_voxels = len(moose.vec('/model/chem/dend/membrane').n)
           Mvec = np.asarray( [0] * num_voxels )
           update_mem = membrane_old.copy()
           old_mem_reinit = membrane_old.copy()
           Mvec[dome_start_voxels[cv] : dome_end_voxels[cv] + 1] = 1
           this_phitot = sum(membrane_old[dome_start_voxels[cv]:dome_end_voxels[cv] + 1])
           print("Mvec: ", Mvec)
           print(dome_start_voxels[cv], dome_end_voxels[cv])
           for i in range(num_voxels):
              update_mem[i] = Mvec[i] * this_phitot / (dome_end_voxels[cv] - dome_start_voxels[cv] + 1)
              old_mem_reinit[i] = (1 - Mvec[i]) * membrane_old[i]
           moose.vec( '/model/chem/dend/membrane' ).n = np.asarray(update_mem) + np.asarray(old_mem_reinit)    
           #moose.vec( '/model/chem/dend/membrane' ).n += np.asarray(update_mem)
        print("After magic diffusion:")
        print("Membrane", moose.vec( '/model/chem/dend/membrane' ).n)
        sum_membrane = sum(moose.vec( '/model/chem/dend/membrane' ).n)
        sum_cytosol = sum(moose.vec( '/model/chem/dend/cytosol' ).n)
        sum_total = sum_membrane + sum_cytosol
        print("Total molecules: ", sum_total)

def delete_spine(delete_voxels):
        membrane_old = moose.vec('/model/chem/dend/membrane').n
        cytosol_old = moose.vec('/model/chem/dend/cytosol').n
        marker = moose.vec('/model/chem/dend/marker').n
        phitot = 0
        for cv in delete_voxels:
            phitot = phitot + membrane_old[cv]
            cytosol_old[cv] = cytosol_old[cv] + membrane_old[cv]
            membrane_old[cv] = 0
            marker[cv] = 0
        moose.vec('/model/chem/dend/membrane').n = membrane_old
        moose.vec('/model/chem/dend/cytosol').n = cytosol_old
        moose.vec('/model/chem/dend/marker').n = marker

def contiguous_chunk(list_item):
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

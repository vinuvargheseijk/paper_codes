import totalEnergy as TE
import numpy as np
import matplotlib.pyplot as plt
import single_data_extract
import aggregation_energy as AE
import pandas as pd

scaling = 1e16
spines = [1]
sum_energy = []
for ns in spines:
   ds = single_data_extract.data_extract(ns)
   phi_list, rm_list, theta_list, Hmean_list, mConc_list, cConc_list, rp_list = ds.phi_list, ds.rm_list, ds.theta_list, ds.Hmean_list, ds.mConc_list, ds.cConc_list, ds.rp_list

   phi_list = np.asarray(phi_list)
   rm_list = np.asarray(rm_list)
   theta_list = np.asarray(theta_list)
   Hmean_list = np.asarray(Hmean_list)
   mConc_list = np.asarray(mConc_list)
   cConc_list = np.asarray(cConc_list)
   rp_list = np.asarray(rp_list)


   area = 2 * np.pi * rm_list**2 * (1 - np.cos(theta_list))
   phitot_list = phi_list * area * TE.phisat 
   print(phitot_list)
   TE.Fminus( rm_list, theta_list, phitot_list, TE.phisat, TE.khat, TE.fs, TE.k )
   entropy_energy = area * TE.entropy(area, phitot_list, TE.phisat, TE.khat, TE.fs, TE.k) * scaling
   aggregation_energy_list = AE.aggregation_energy(rm_list, theta_list, phitot_list, TE.phisat, TE.k_agg) * scaling
   mismatch_list = TE.mismatch( rm_list, phitot_list, TE.phisat, TE.khat, TE.fs, TE.k ) * scaling
   bending_dome = ( area * TE.k / ( 2 * rm_list**2 ) ) * scaling
   bending_saddle = TE.Fplus(rp_list, rm_list, theta_list, TE.khat, TE.fs, TE.k) * scaling
   tension_energy = TE.membrane_tension(rm_list, theta_list, rp_list, phitot_list, TE.phisat, 1, TE.khat, TE.fs, TE.k, 0, 0) * scaling
   mechanical_potential = AE.mech_potential( Hmean_list, phi_list, TE. k_agg) * scaling
   #potential_energy_full = phitot_list * TE.full_potential(phitot_list, rm_list, theta_list, cConc_list, mConc_list) * scaling
   potential_energy_chem = phitot_list * TE.chem_pot_partition(phitot_list, cConc_list, rm_list * 1e6 , theta_list) * scaling
   disk_energy = TE.Fzero( rp_list, rm_list, theta_list, 1, TE.khat, TE.fs, TE.k, 0, 0) * scaling
   constraint_energy = TE.area_constraint(rp_list , rm_list, theta_list, TE.dendDia, TE.fs, phitot_list, area) * scaling

   sum_energies = entropy_energy[-1] + aggregation_energy_list[-1] + mismatch_list[-1] + bending_dome[-1] + bending_saddle[-1] + tension_energy[-1] + potential_energy_chem[-1] + constraint_energy[-1]
   print("Sum of all: ", sum_energies)
   
   sum_energies_abs = np.abs(entropy_energy[-1]) + np.abs(aggregation_energy_list[-1]) + np.abs(mismatch_list[-1]) + np.abs(bending_dome[-1]) + np.abs(bending_saddle[-1]) + np.abs(tension_energy[-1]) + np.abs(potential_energy_chem[-1]) + np.abs(constraint_energy[-1])
   
   df_energies = pd.DataFrame(columns = ["Ee", "Ea", "Em", "Eb", "Et", "Ep", "Ec"])

   Ee_frac = np.abs(entropy_energy[-1]) / sum_energies_abs
   df_energies["Ee"] = [Ee_frac]
   Ea_frac = np.abs(aggregation_energy_list[-1]) / sum_energies_abs
   df_energies["Ea"] = [Ea_frac]
   Em_frac = np.abs(mismatch_list[-1]) / sum_energies_abs
   df_energies["Em"] = [Em_frac]
   Ebd_frac = np.abs(bending_dome[-1]) / sum_energies_abs
   Ebs_frac = np.abs(bending_saddle[-1]) / sum_energies_abs
   df_energies["Eb"] = [Ebs_frac + Ebd_frac]
   Et_frac = np.abs(tension_energy[-1]) / sum_energies_abs
   df_energies["Et"] = [Et_frac]
   Ep_frac = np.abs(potential_energy_chem[-1]) / sum_energies_abs
   df_energies["Ep"] = [Ep_frac]
   Ec_frac = np.abs(constraint_energy[-1]) / sum_energies_abs
   df_energies["Ec"] = [Ec_frac]
    
   sum_frac_data = 0
   for i in df_energies.iloc[0]:
       sum_frac_data = sum_frac_data + float(i)
   print("SUM FRAC: ", sum_frac_data)    
   df_energies["sum"] = [sum_frac_data]
   print(df_energies)
   df_energies.to_csv("./energy_fraction.csv")

   sum_frac = Ee_frac + Ea_frac + Em_frac + Ebd_frac + Ebs_frac + Et_frac + Ep_frac + Ec_frac



   print("Sum fractions: ", sum_frac)

   plt.plot(mismatch_list, label = "mismatch")
   plt.plot(entropy_energy, label = "Entropy")
   plt.plot(aggregation_energy_list, label = "Aggregation_energy")
   plt.plot(bending_saddle + bending_dome, label = "Bending_total")
   plt.plot(tension_energy, label = "Tension")
   plt.plot(mechanical_potential, label = "Mech. Potential")
   plt.plot(potential_energy_chem, label = "Chem. Potential")
   plt.plot(disk_energy, label = "disk energy")
   plt.plot(constraint_energy, label = "Constraint")
   plt.legend()
   plt.figure(1001)
   totalEnergy = TE.total_energy([ rm_list * 1e6, theta_list, rp_list * 1e6 ], phitot_list, TE.phisat, TE.khat, TE.fs, TE.k, 1, 0, 0, 0, cConc_list)
   plt.title("Total Energy" + " spine: " + str(ns))
   plt.plot(totalEnergy)
   plt.legend()
   plt.show()
   sum_energy.append( list(totalEnergy) )
plt.figure(1002)   
max_length = 0
min_length = 1e8
for ns in range(len(sum_energy)):
    if len( sum_energy[ns] ) > max_length:
        max_length = len( sum_energy[ns] )
    if len( sum_energy[ns] ) < min_length:
        min_length = len( sum_energy[ns] )


#total = np.asarray( [0] * max_length )
for ns in range(len(sum_energy)):
    diff_length = max_length - len(sum_energy[ns])
    temp_list = [0] * diff_length
    sum_energy[ns].extend( temp_list )
    plt.plot( sum_energy[ns][0:min_length], label = str( spines[ns] ) )
    #total = total + np.asarray( sum_energy[ns] )
#plt.plot(total, label = 'Sum') 
plt.legend()
plt.show()






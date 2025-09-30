import pandas as pd
import numpy as np
import dendShape
import writeXML
import matplotlib.pyplot as plt
import totalEnergy
#import chemModel as CM
import readXML as RX
import xml.etree.ElementTree as ET
import totalEnergy as TE
import tag
import seaborn as sns
from scipy import optimize

sns.set_context("poster")


mu0 = totalEnergy.calc_mu0()
RT = totalEnergy.RT
khat = totalEnergy.khat
KBT = totalEnergy.KBT
Cp = totalEnergy.Cp
phisat = totalEnergy.phisat
k_agg = totalEnergy.k_agg
Kb = tag.Kb
part_coeff = totalEnergy.part_coeff
diffL = totalEnergy.diffL
shape_dt = 0.1
c_dt = shape_dt
plot_interval = 1
diffL = totalEnergy.diffL
Length = tag.Length * 1e-6
num_voxels = int(Length / diffL)
V_voxel = TE.V_cyt / num_voxels
Rc = 0.5 * totalEnergy.dendDia 
fs = totalEnergy.fs
print(tag.tag_string)
def remove_nan(list_v):
    nonNan = []
    for i in list_v:
         if np.isnan(i) == False:
            nonNan.append(i)
    return nonNan

def all_locations():
    sstart = pd.read_csv(tag.tag_string + "sstart.csv")
    hmean = pd.read_csv(tag.tag_string + "hmean.csv")
    theta = pd.read_csv(tag.tag_string + "theta.csv")
    rp = pd.read_csv(tag.tag_string + "rp.csv")
    max_locations = []
    max_spines = 0
    for i in range(len(sstart)):
        if len(remove_nan(sstart.iloc[i][1:])) > max_spines:
           max_spines = len(remove_nan(list(sstart.iloc[i][1:])))
           max_locations = remove_nan(list(sstart.iloc[i][1:]))
           index = i
    hmean_max = remove_nan(hmean.iloc[index][1:])
    theta_max = remove_nan(theta.iloc[index][1:])
    rp_max = remove_nan(rp.iloc[index][1:])
    dome_start = np.asarray(max_locations) + np.asarray(rp_max) * np.sin(np.asarray(theta_max))
    dome_end = dome_start + 2 * np.asarray(1 / np.asarray(hmean_max)) * np.sin(np.asarray(theta_max))
    return index, max_locations, dome_start, dome_end

index, max_locations, dome_start, dome_end = all_locations()
df_locations = pd.DataFrame()
df_locations["locations"] = max_locations
df_locations["dome_start"] = dome_start
df_locations["dome_end"] = dome_end
df_locations.to_csv(tag.tag_string + "MaxLocations.csv")
print("INDEX, locations: ", index, max_locations)
print("Dome: ", dome_start)

#df_exp = pd.read_csv("./Ca_data_singleAp_Chen_2013.csv")
def calc_conc_voxel(numbers):
    return 1e3 * np.asarray(numbers) / (TE.V_voxel * TE.Na)

def normalize_exp_sim(data):
    return np.asarray(data) / max(data)

def plot_scale( ax, data ):
    if True in np.isnan(data):
        print("nan found")
    if len(data) == 0:
        print("")
    else:    
        ticks = np.arange(0, max(data), 250)  
        print(ticks)
        ax.set_yticks( ticks )
    
def pull_cost(rs, rm, theta, alpha, frac):
    rs = rs * 1e-6
    rm = rm * 1e-6
    yr = rs * (1 - np.cos(alpha))
    print(rs, rm, yr)
    area_saddle = frac * totalEnergy.calc_saddle_area(rm, theta, rs)
    print(area_saddle)
    frac_area_energy = area_saddle * ( (yr - Rc)**2 / Rc**2 ) * fs
    print("FAE: ", frac_area_energy)
    return frac_area_energy / KBT

def tension_relax(rs, rm, theta, alpha, frac):
    rs = rs * 1e-6
    rm = rm * 1e-6
    frac_area = frac * totalEnergy.calc_saddle_area(rm, theta, rs)
    return frac_area * fs / KBT

def sum_N2(rm, theta, rp, phitot, conc, frac):
    total_energy = 0
    frac_gain = 0
    for n in range(len(rm)):
        total_energy = total_energy + totalEnergy.total_energy([rm[n] * theta[n], theta[n], rp[n]], phitot[n], conc[n]) 
        frac_gain = frac_gain + totalEnergy.total_energy([rm[n] * theta[n], theta[n], rp[n]], phitot[n], conc[n]) - tension_relax(rp[n], rm[n], theta[n], theta[n], frac) + pull_cost(rp[n], rm[n], theta[n], theta[n], frac)
        print(rp[n], rm[n], theta[n])

        print(total_energy)
        print(frac_gain)
    return total_energy, frac_energy    


class data_extract:

  def plot_ref(iterations, ax1):
      num_refs = len(converged_phitots)
      for nr in range(num_refs):
         plot_phitot = [converged_phitots["final_phitot"][nr]] * len(globals()["phitot" + str(1) + "_all"]) 
         ax1.plot(iterations, plot_phitot, label = str(converged_phitots["phi_entire"][nr]))
      ax1.legend()


  def data_save_lists():
    for ns in range(30):
      globals()["phi" + str(ns+1) + "_all"] = []
      globals()["phitot" + str(ns+1) + "_all"] = []
      globals()["sum_phitot_all"] = []
      globals()["theta" + str(ns+1) + "_all"] = []
      globals()["rm" + str(ns+1) + "_all"] = []
      globals()["hmean" + str(ns+1) + "_all"] = []
      globals()["rp" + str(ns+1) + "_all"] = []
      globals()["mConc" + str(ns+1) + "_all"] = []
      globals()["cConc" + str(ns+1) + "_all"] = []
      globals()["voxel_count" + str(ns) + "_all"] = []
      globals()["energy" + str(ns) + "_all"] = []
      globals()["Ae" + str(ns) + "_all"] = []
      globals()["Ee" + str(ns) + "_all"] = []
      globals()["Me" + str(ns) + "_all"] = []
      globals()["Be" + str(ns) + "_all"] = []
      globals()["Ce" + str(ns) + "_all"] = []
    tag_string = tag.tag_string
    phi = pd.read_csv(tag_string + "phi.csv")
    theta = pd.read_csv(tag_string + "theta.csv")
    saddle_start = pd.read_csv(tag_string + "sstart.csv")
    hmean = pd.read_csv(tag_string + "hmean.csv")
    rp = pd.read_csv(tag_string + "rp.csv")
    mConc = pd.read_csv(tag_string + "mConc.csv")
    cConc = pd.read_csv(tag_string + "cConc.csv")




    y_save = []
    x_save = []
    tEn_list = []
    for i in range(len(phi)):
      theta_list = []
      Hmean_list = []
      Rp_list = []
      saddle_start_list = []
      loc_index = 0
      sum_phitot = 0
      spine_wise_sum_energy = 0
      for j in range(1, len(phi.iloc[i])):
        if np.isnan(theta.iloc[i][j]) == True:
          globals()["Ae" + str(j) + "_all"].append( 0 )
          globals()["Ee" + str(j) + "_all"].append( 0 )
          globals()["Me" + str(j) + "_all"].append( 0 )
          globals()["Be" + str(j) + "_all"].append( 0 )
          globals()["Ce" + str(j) + "_all"].append( 0 )
            

        if np.isnan(theta.iloc[i][j]) == False:
          theta_list.append( theta.iloc[i][j] )
          Hmean_list.append( hmean.iloc[i][j] )
          Rp_list.append( rp.iloc[i][j] )
          saddle_start_list.append( saddle_start.iloc[i][j] )
          loc_index = loc_index + 1
          globals()["phi" + str(j) + "_all"].append(phi.iloc[i][j])
          area = 2 * np.pi * ( 1 / hmean.iloc[i][j] )**2 * (1 - np.cos(theta.iloc[i][j]))
          phitot = phi.iloc[i][j] * area * phisat
          #current_energy = optimize.minimize( TE.total_energy, x0=TE.x0, bounds = (TE.rmbounds,TE.tbounds,TE.rpbounds),args=(phitot,TE.phisat, TE.khat, TE.fs, TE.k, 1, 0, 0, 0, cConc.iloc[i][j]),method = 'L-BFGS-B', options = {'ftol':1e-12} )
          alength_i = (1 / hmean.iloc[i][j]) * theta.iloc[i][j] * 1e6
          current_energy = TE.total_energy( [alength_i, theta.iloc[i][j], rp.iloc[i][j] * 1e6], phitot, cConc.iloc[i][j] )
          spine_wise_sum_energy = spine_wise_sum_energy + current_energy

          sum_phitot = sum_phitot + phitot
          globals()["phitot" + str(j) + "_all"].append(phitot)

          globals()["rm" + str(j) + "_all"].append( 1 / hmean.iloc[i][j] )
          globals()["rp" + str(j) + "_all"].append( rp.iloc[i][j] )
          globals()["theta" + str(j) + "_all"].append(theta.iloc[i][j])
          globals()["hmean" + str(j) + "_all"].append(hmean.iloc[i][j])
          globals()["mConc" + str(j) + "_all"].append(mConc.iloc[i][j])
          globals()["cConc" + str(j) + "_all"].append(cConc.iloc[i][j])
          globals()["energy" + str(j) + "_all"].append(current_energy)
          globals()["Ae" + str(j) + "_all"].append(TE.aggregation_energy(phi.iloc[i][j], phitot, area) / TE.KBT )
          globals()["Ee" + str(j) + "_all"].append(TE.entropy(phi.iloc[i][j], phitot, area, TE.khat ) / TE.KBT )
          globals()["Me" + str(j) + "_all"].append(TE.mismatch( 1/hmean.iloc[i][j], phi.iloc[i][j], phitot, area ) / TE.KBT)
          globals()["Be" + str(j) + "_all"].append(TE.dome_bending( 1/hmean.iloc[i][j], phi.iloc[i][j], area) / TE.KBT)
          globals()["Ce" + str(j) + "_all"].append(TE.chem_pot_partition(phitot, cConc.iloc[i][j], 1/hmean.iloc[i][j], mu0, area)/ TE.KBT)
      globals()["sum_phitot_all"].append(sum_phitot)
      tEn_list.append(spine_wise_sum_energy)
      if len(theta_list) > 0:    
         x, y = dendShape.get_shape(saddle_start_list, theta_list, Hmean_list, Rp_list)
         x_save.append(x)
         y_save.append(y)

    try:
      rm_last = np.asarray([globals()["rm" + str(1) + "_all"][-1], globals()["rm" + str(2) + "_all"][-1]]) * 1e6
      theta_last = [ globals()["theta" + str(1) + "_all"][-1],  globals()["theta" + str(2) + "_all"][-1] ]
      rp_last = np.asarray( [ globals()["rp" + str(1) + "_all"][-1],  globals()["rp" + str(2) + "_all"][-1] ] ) * 1e6
      phitot_last = [ globals()["phitot" + str(1) + "_all"][-1],  globals()["phitot" + str(2) + "_all"][-1] ]
      conc_last = [ globals()["cConc" + str(1) + "_all"][-1],  globals()["cConc" + str(2) + "_all"][-1] ]
      print(rm_last, theta_last, rp_last, phitot_last, conc_last)
      summ_of_F2, merge_cost  = sum_N2(rm_last, theta_last, rp_last, phitot_last, conc_last, 0.5)
      print("Sum of first two spines: ", summ_of_F2, merge_cost)                    
    except Exception as excpt:
      print(type(excpt).__name__)  
      print("N<2")                    

    for ne in range(1, 3):
        df_e = pd.DataFrame()
        try:
           print(ne) 
           df_e["Ae"+str(ne)] = globals()["Ae" + str(ne) + "_all"]
           df_e["Ee"+str(ne)] = globals()["Ee" + str(ne) + "_all"]
           df_e["Me"+str(ne)] = globals()["Me" + str(ne) + "_all"]
           df_e["Be"+str(ne)] = globals()["Be" + str(ne) + "_all"]
           df_e["Ce"+str(ne)] = globals()["Ce" + str(ne) + "_all"]
           print("N = " + str(ne) + " inserted")
        except Exception as e:
           print(e) 
           print("Number of spines exceeded")
        print(df_e)       
        df_e.to_csv("./" + tag.tag_string + "energies" + str(ne) + ".csv")
    t_axis = pd.read_csv(tag_string + 'chemTimes.csv')["time"] 
    x_axis = pd.read_csv(tag_string + 'shapeTimes.csv')["0"]
    writeXML.writeXML_td(x_save,tag_string + "X.xml")
    writeXML.writeXML_td(y_save,tag_string + "Y.xml")
    fig, ( (ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8), (ax9, ax10) ) = plt.subplots(5, 2)
    ax21 = ax2.twinx()
    ax1.set_title("$\phi_{tot}$ solid, cyt Conc. dashed")
    ax11 = ax1.twinx()
    ax61 = ax6.twinx()
    ax11.set_ylabel( "Cytosol Conc." )
    ax1.set_ylabel( "$\phi_{tot}$" )
    ax2.set_title("$R^{d}$")
    ax9.set_title("merge cost")
    if voxel_detect == False:
        ax9.set_title("$\Theta$")
    ax4.set_title("Mechanical potential")
    ax5.set_title("Chemical potential")
    ax6.set_title("Rate Kf (solid), Residue (dashed)")
    ax7.set_title("Total Energy")
    ax71 = ax7.twinx()
    #ax71.set_ylabel( "Av. [Cyt]" )
    if plot_ref_flag:
       iterations = np.arange( 0, len(globals()["phitot" + str(1) + "_all"]) )
       data_extract.plot_ref(iterations, ax1)
    #ax7.plot( globals()["sum_phitot_all"], label = "Total $\phi_{tot}$" )
    #plot_scale( ax7, globals()["sum_phitot_all"])
    ax7.plot(x_axis, tEn_list)
    #ax7.set_ylim(0, 100)
    for ns in range(1, 10):
     len_data = len(globals()["phitot" + str(ns) + "_all"])
     x_axis = x_axis[0:len_data]
     ax1.plot( x_axis, globals()["phitot" + str(ns) + "_all"] )
     #plot_scale( ax1, globals()["phitot" + str(ns) + "_all"] )
     ax11.plot( x_axis, globals()["cConc" + str(ns) + "_all"], linestyle = "-." )
     ax2.plot( x_axis, globals()["rm" + str(ns) + "_all"], 'r')
     ax2.set_ylabel("$R^{d}$")
     ax21.set_ylabel("$R^{s}$")
     ax21.plot( x_axis, globals()["rp" + str(ns) + "_all"], 'g')
     #ax7.plot( globals()["energy" + str(j) + "_all"], label = "Energy" )
     ax8.plot( x_axis, globals()["phi" + str(ns) + "_all"])
     ax8.set_title("$\phi$")
     rp_array = np.asarray( globals()["rp" + str(ns) + "_all"] )
     sum_rad = rp_array + np.asarray( globals()["rm" + str(ns) + "_all"] )
     r0_array = sum_rad * np.sin( np.asarray( globals()["theta" + str(ns) + "_all"] ) )
     ax2.plot( x_axis, r0_array, linestyle = '-.' ) 
     if voxel_detect == False:
        ax3.plot(x_axis, globals()["theta" + str(ns) + "_all"] )
     #ax4.plot( globals()["hmean" + str(ns) + "_all"] )
     hmean_array = np.asarray( globals()["hmean" + str(ns) + "_all"] )
     rm_array = np.asarray( globals()["rm" + str(ns) + "_all"] )
     theta_array = np.asarray( globals()["theta" + str(ns) + "_all"] )
     phi_array = np.asarray( globals()["phi" + str(ns) + "_all"] )
     concRatio = np.asarray(globals()["mConc" + str(ns) + "_all"]) / np.asarray(globals()["cConc" + str(ns) + "_all"])
     conc_array = np.asarray(globals()["cConc" + str(ns) + "_all"])
     mconc_array = np.asarray(globals()["mConc" + str(ns) + "_all"])
     phitot_array = np.asarray(globals()["phitot" + str(ns) + "_all"])
     logConcRatio = np.log( concRatio )
     secondTerm = RT * logConcRatio
     firstTerm = mu0
     #chemicalPotential = firstTerm + secondTerm
     areaArray = 2 * np.pi * rm_array**2 * (1 - np.cos(theta_array))
     chemicalPotential = TE.chem_pot_partition(phitot_array, conc_array, rm_array * 1e6, firstTerm, areaArray) / phitot_array
     ax5.plot( x_axis, chemicalPotential / RT )
     mechanical_potential = TE.mech_potential(hmean_array, phi_array)
     ax4.plot( x_axis, mechanical_potential / RT )
     kf_values = Kb * part_coeff * np.exp( - ( mechanical_potential + chemicalPotential ) / RT )
     ax6.plot( x_axis, kf_values )
     residue_steady = kf_values * conc_array - Kb *  mconc_array
     ax61.plot( x_axis, residue_steady, linestyle = '-.' )


    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.7, hspace=0.85)
    if voxel_detect == True:
       voxel_count = [] 
       filename = tag_string + "marker.xml"
       filename_m = tag_string + "membrane.xml"
       tree1 = ET.parse(filename)
       tree = [tree1]
       sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
       tree2 = ET.parse(filename_m)
       tree = [tree2]
       sum_an_m,max_an_y_m,an_y_m=RX.plotXML(filename_m,tree)
       sum_mem_marker = []
       for l in range(0, len(an_y)):
          vox_temp_count = 0
          ns = 1
          nz = np.nonzero( np.asarray(an_y[l]) )[0]
          temp = 0
          for z in range( len(nz) - 1 ):
              if nz[z] + 1 == nz[z + 1]:
                  vox_temp_count = vox_temp_count + 1
                  temp = temp + an_y_m[l][nz[z]]
              else:
                  globals()["voxel_count" + str(ns) + "_all"].append( vox_temp_count + 1 )
                  temp = temp + an_y_m[l][nz[z] + 1]
                  vox_temp_count = 0
                  ns = ns + 1
          sum_mem_marker.append(temp)     
       #ax1.plot(sum_mem_marker, label = 'under dome') 
       #ax1.legend()
       for ns in range(1, 10):
           lv = len(globals()["voxel_count" + str(ns) + "_all"])
           range_iter = np.arange(len(an_y) - lv, len(an_y), 1) 
           ax3.plot( range_iter, globals()["voxel_count" + str(ns) + "_all"] )
    if conserv_check == True:   
       membraneXML = tag_string + "membrane.xml"
       cytosolXML = tag_string + "cytosol.xml"
       check_files = [membraneXML, cytosolXML]
       sum_mols = [0] * len(an_y)
       for c in check_files:
         filename = c
         tree1 = ET.parse(filename)
         tree = [tree1]
         sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
         mol_count = []
         print(num_voxels, len(an_y[-1]))
         for l in range(len(an_y)):
             mol_count.append( sum( an_y[l] ) )
         label = c.replace( tag_string, "" )    
         label = label.replace( ".xml", "" )    
         """
         ax8.plot( mol_count, label = label )
         
         if c == cytosolXML:
             ax71.plot( np.asarray( mol_count ) / ( TE.Na * TE.V_cyt ), 'r', label = "[Cyt]" )
         """    
         sum_mols = np.asarray( sum_mols ) + np.asarray( mol_count )
       CaXML = tag_string + "Ca.xml"
       RacXML = tag_string + "Rac.xml"
       thrXML = tag_string + "thr286.xml"
       check_files = [CaXML]
       names = ["Ca"]
       mol_name = 0
       ax11 = ax10.twinx()
       #norm_exp = normalize_exp_sim(df_exp["ca"])
       #time_exp = df_exp["time"]
       #ax10.plot(time_exp, norm_exp)
       for c in check_files:
         filename = c
         tree1 = ET.parse(filename)
         tree = [tree1]
         sum_an,max_an_y,an_y=RX.plotXML(filename,tree)
         if names[mol_name] == "Rac":
            ax11.plot( t_axis, calc_conc_voxel(max_an_y), 'r', label = names[mol_name] )
         if names[mol_name] == "Ca":
            ax10.plot( t_axis, calc_conc_voxel(max_an_y), 'b', label = names[mol_name] )
         if names[mol_name] == "thr286":
            ax9.plot( t_axis, calc_conc_voxel(max_an_y), 'b', label = names[mol_name] )
         mol_name = mol_name + 1
         ax10.set_title( "Rac/Ca") 
         ax10.legend()
         ax11.legend()
       """
       ax71.legend()
       """
       ax7.legend()
    
    plt.show()

plot_ref_flag = False
voxel_detect = True
conserv_check = True

if plot_ref_flag:
   converged_phitots = pd.read_csv("converged_values.csv")
       
   
data_extract.data_save_lists()


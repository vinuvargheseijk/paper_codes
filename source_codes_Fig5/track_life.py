import numpy as np
import matplotlib.pyplot as plt
import statistics
import pandas as pd

figH = plt.figure(figsize = (10, 8))
ax1H = plt.subplot2grid((2, 3), (0, 0), colspan = 1) 
ax11H = plt.subplot2grid((2, 3), (0, 1), colspan = 1) 
ax12H = plt.subplot2grid((2, 3), (0, 2), colspan = 1) 
ax21H = plt.subplot2grid((2, 3), (1, 0), colspan = 1) 
ax23H = plt.subplot2grid((2, 3), (1, 2), colspan = 1) 
ax11Ht = ax11H.twinx()
ax12Ht = ax12H.twinx()
ax21Ht = ax21H.twinx()
trials = list(np.arange(0, 20, 1))
start_spine = 11e-6
end_spine = 39e-6
start_sync_range = 24e-6
NstimL = 10.0
spacing = 1.0
shape_dt = 0.2
bfreq = 20.0
end_sync_range = start_sync_range + NstimL * spacing * 1e-6 + 1e-6



def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

def spot_boundary_spine(coordinate):
    proximity = False
    if abs(coordinate - start_spine) < 1e-6:
       proximity = True
    if abs(coordinate - end_spine) < 1e-6:
       proximity = True
    return proximity

def get_sync_spines(coords):
    ss = coords[np.where(coords > start_sync_range)]
    ss1 = ss[np.where(ss < end_sync_range)]
    return ss1

def get_out_sync_spines(coords):
    iss = list(coords[np.where(coords < start_sync_range)])
    for c in iss:
        proximal = spot_boundary_spine(c)
        if proximal:
           iss.remove(c)
    iss1 = list(coords[np.where(coords > end_sync_range)])
    for c in iss1:
        proximal = spot_boundary_spine(c)
        if proximal:
           iss1.remove(c)
    return iss, iss1

def isd_calc(saddle_coordinates):
    if len(saddle_coordinates) > 1:
       isd = (saddle_coordinates[-1] - saddle_coordinates[0]) / (len(saddle_coordinates) - 1)
    else:
       isd = 0 
    return isd

def spine_density(saddle_coordinates, region):
    if len(saddle_coordinates) > 0:
       if region == "Rs":
          spine_density = len(saddle_coordinates) / 10
       if region == "Rns":
          spine_density = len(saddle_coordinates) / 20
    else:
       spine_density = 0
    return spine_density

def order_list(lst):
    for i in range(len(lst)):
       for j in range(len(lst)):
          if lst[i][0] < lst[j][0]:
             temp = lst[i]
             lst[i] = lst[j]
             lst[j] = temp
    return lst

def get_bins(coords, coord_bins, time):
    for cd in coords:
       if len(coord_bins) > 0:
          bin_add_flag = True
          for cbins in coord_bins:
              if cd >= cbins[0] and cd <= cbins[1]:
                 cd_index = coord_bins.index(cbins)
                 bin_add_flag = False
          if bin_add_flag:
             coord_bins.append([cd - 0.5e-6, cd + 0.5e-6])
       else:         
          coord_bins.append([cd - 0.5e-6, cd + 0.5e-6])
    return order_list(coord_bins)

def avg_lifetime(dict_spines):
    spines_index = []
    life_index = []
    for sp in range(len(dict_spines)):
       if dict_spines[sp] > 0:
         spines_index.append(sp)
         life_index.append(dict_spines[sp])
    print("LIFE INDEX: ", life_index)
    if len(life_index) > 0:
        mean_life = statistics.mean(life_index)
    else:
        mean_life = 0
    return mean_life

def get_life(enz_kcat):
  stim_patterns = ["stim"]
  parString = "NstimL_10.0_kcat_enz_tiam_" + str(enz_kcat) + "_spine_spacing_30.0_spacing_1.0_bfreq_" + str(bfreq) + "_coactive_1.0_mixed_True_trial_0"
  sync_isd_list = []  #Mean along time
  out_sync_isd_list = []
  sync_life_time = []
  out_sync_life_time = []
  sync_sp_density_list = []
  out_sync_sp_density_list = []
  missing_trials = []
  for tnum in range(len(trials)):
    try: 
       life_dict = {}
       life_dict_sync = {}
       life_dict_out_sync = {}
       print(parString)
       print(tnum)
       phi = pd.read_csv("./" + parString + "phi.csv")
       sstart = pd.read_csv("./" + parString + "sstart.csv")
       shapeTimes = pd.read_csv("./" + parString + "shapeTimes.csv")
       coord_bins = []
       df = pd.DataFrame()
       ISD_sync_spines = []
       ISD_out_sync_spines = []
       sp_density_sync_spines = []
       sp_density_out_sync_spines = []
       for it in range(len(shapeTimes)):
          coords = remove_nan(sstart.iloc[it][2:-1])
          coord_bins = get_bins(coords, coord_bins, it)
          #print(coords, coord_bins) 
          sync_coords = get_sync_spines(np.asarray(coords))
          temp_sync_isd = isd_calc(sync_coords)
          ISD_sync_spines.append(temp_sync_isd)
          out_sync_coords_left, out_sync_coords_right = get_out_sync_spines(np.asarray(coords))    
          ISD_out_sync_spines.append(isd_calc(out_sync_coords_left) + isd_calc(out_sync_coords_right))
          sp_density_sync_spines.append(spine_density(sync_coords, "Rs"))
          sp_density_out_sync_spines.append(spine_density(out_sync_coords_left + out_sync_coords_right, "Rns"))
       
       #final_cords = np.asarray(remove_nan( list(sstart.iloc[-1][1:]) ))
       #sync_spine_cords = get_sync_spines(final_cords)
       #out_sync_spines = get_out_sync_spines(final_cords)
       sync_isd_list.append( statistics.mean(ISD_sync_spines) )
       out_sync_isd_list.append( statistics.mean(ISD_out_sync_spines) )
       sync_sp_density_list.append( statistics.mean(sp_density_sync_spines) )
       out_sync_sp_density_list.append( statistics.mean(sp_density_out_sync_spines) )

       for msp in range(30):
          life_dict[msp] = 0
          life_dict_sync[msp] = 0
          life_dict_out_sync[msp] = 0
       for it in range(1, len(shapeTimes)):
          coords = remove_nan(sstart.iloc[it][1:])
          for cd in coords:
             for cbins in coord_bins:
                if cd >= cbins[0] and cd <= cbins[1]:
                  bin_index = coord_bins.index(cbins)
                  life_dict[bin_index] = life_dict[bin_index] + (shapeTimes["0"][it] - shapeTimes["0"][it - 1])
                  if cd >= start_sync_range and cd <= end_sync_range:
                      life_dict_sync[bin_index] = life_dict_sync[bin_index] + (shapeTimes["0"][it] - shapeTimes["0"][it - 1])
                  if cd < start_sync_range - 1e-6 or cd > end_sync_range + 1e-6:
                     if cd > 10.5e-6 and cd < 39.5e-6:
                         life_dict_out_sync[bin_index] = life_dict_out_sync[bin_index] + (shapeTimes["0"][it] - shapeTimes["0"][it - 1])
    except FileNotFoundError:
          print("kcat: " + str(enz_kcat) + ", trial: " + str(tnum) + " not found")
          missing_trials.append(tnum)
    if tnum < len(trials) - 1:
      parString = parString.replace("trial_" + str(trials[tnum]), "trial_" + str(trials[tnum + 1]))
    avg_life_sync = avg_lifetime(life_dict_sync)
    avg_life_out_sync = avg_lifetime(life_dict_out_sync)
    sync_life_time.append(avg_life_sync)
    out_sync_life_time.append(avg_life_out_sync)
  return sync_life_time, out_sync_life_time, sync_isd_list, out_sync_isd_list, sync_sp_density_list, out_sync_sp_density_list

def get_size(enz_kcat):
  parString = "NstimL_10.0_kcat_enz_tiam_" + str(enz_kcat) + "_spine_spacing_30.0_spacing_1.0_bfreq_" + str(bfreq) + "_coactive_1.0_mixed_True_trial_0"
  avg_sync_size = []
  avg_out_of_sync_size = []
  avg_area_sync = []
  avg_area_out_of_sync = []
  missing_trials = []
  for tnum in range(len(trials)):
     try:
        theta = pd.read_csv("./" + parString + "theta.csv")
        hmean = pd.read_csv("./" + parString + "hmean.csv")
        rp = pd.read_csv("./" + parString + "rp.csv")
        sstart = pd.read_csv("./" + parString + "sstart.csv")
        final_hmean = remove_nan(hmean.iloc[-1][2:-1])
        final_theta = remove_nan(theta.iloc[-1][2:-1])
        final_rp = remove_nan(rp.iloc[-1][2:-1])
        coords = remove_nan(sstart.iloc[-1][2:-1])
        size_sync = []
        size_out_of_sync = []
        area_sync = []
        area_out_of_sync = []
        for cd in range(len(coords)):
          if coords[cd] >= start_sync_range and coords[cd] <= end_sync_range:
             print("Coords in SYNC: ", coords[cd])
             temp_size =  2 * (1 / final_hmean[cd] + final_rp[cd]) * np.sin(final_theta[cd])
             temp_area = np.pi * (1 / final_hmean[cd])**2 * (1 - np.cos(final_theta[cd]))
             size_sync.append( temp_size )
             area_sync.append( temp_area )
          if coords[cd] < start_sync_range or coords[cd] > end_sync_range:
             if coords[cd] < end_spine:
                print("Coords in ASYNC: ", coords[cd])
                size_out_of_sync.append( 2 * (1 / final_hmean[cd] + final_rp[cd]) * np.sin(final_theta[cd]) )
                area_out_of_sync.append( np.pi * (1 / final_hmean[cd])**2 * (1 - np.cos(final_theta[cd])) )
        print("AREA SYNC: ", area_sync)
        print("AREA OutSYNC: ", area_out_of_sync)
        if len(size_sync) > 0:
           avg_sync_size.append(np.average(size_sync))
           avg_area_sync.append(np.average(area_sync))
        else:
          avg_sync_size.append( 0 )
          avg_area_sync.append( 0 )
        if np.isnan(np.average(area_out_of_sync)):
           avg_out_of_sync_size.append(0)
           avg_area_out_of_sync.append(0)
        else:
           avg_out_of_sync_size.append(np.average(size_out_of_sync))
           avg_area_out_of_sync.append(np.average(area_out_of_sync))
     except FileNotFoundError:
          print("kcat: " + str(enz_kcat) + ", trial: " + str(tnum) + " not found")
          missing_trials.append(tnum)
     if tnum < len(trials) - 1:
        parString = parString.replace("trial_" + str(trials[tnum]), "trial_" + str(trials[tnum + 1]))
  print("KCAT and area: ", enz_kcat, avg_area_sync, avg_area_out_of_sync)
  df_missing = pd.DataFrame()
  df_missing["tmiss"] = missing_trials
  df_missing.to_csv("missingKcatTrial_" + str(enz_kcat) + "_" + str(bfreq) + ".csv")
  return avg_sync_size, avg_out_of_sync_size, avg_area_sync, avg_area_out_of_sync
    
get_size(20.0)

#kcat_list = [20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 200.0]
kcat_list = [20.0, 50.0, 200.0]
for kcat in kcat_list:
    sync_life_time, out_sync_life_time, sync_isd_list, out_sync_isd_list, sync_sp_density_list, out_sync_sp_density_list = get_life(kcat) 
    avg_sync_size, avg_out_of_sync_size, avg_area_sync, avg_area_out_of_sync = get_size(kcat)
    df = pd.DataFrame()
    df["sync_life"] = sync_life_time
    df["out_sync_life"] = out_sync_life_time
    df["sync_isd"] = sync_isd_list
    df["out_sync_isd"] = out_sync_isd_list
    df["sync_size"] = avg_sync_size
    df["out_sync_size"] = avg_out_of_sync_size
    df["sync_area"] = avg_area_sync
    df["out_sync_area"] = avg_area_out_of_sync
    df["sync_sp_density"] = sync_sp_density_list
    df["out_sync_sp_density"] = out_sync_sp_density_list
    df.to_csv("./" + str(kcat) + "_" + str(bfreq) + ".csv")

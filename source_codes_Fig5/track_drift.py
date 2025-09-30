import numpy as np
import matplotlib.pyplot as plt

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
trials = np.arange(0, 20, 1)
start_spine = 10e-6
end_spine = 40e-6
start_sync_range = 24e-6
NstimL = 10.0
spacing = 1.0
end_sync_range = start_sync_range + NstimL * spacing * 1e-6 + 1e-6



def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

def MSD(space, time):
    duration = len(space)
    msd = []
    hist = ax23H.hist(time)
    time_chunks = []
    for i in range(len(hist[0])):
         if hist[0][i] != 0:
           time_chunks.append([hist[1][i], hist[1][i + 1]])
    """
    for dr in range(duration - 1):
        msd.append( ( 1e-3 * (space[dr + 1] - space[dr])**2 / (time[dr + 1] - time[dr]) ) )  #in um^2/ms
        #msd.append( ( 1e-3 * (space[dr + 1] - space[dr]) / (time[dr + 1] - time[dr]) ) )  #in um^2/ms
    """
    for tc in time_chunks:
        for dr in range(duration - 1):
           if time[dr + 1] <= tc[1] and time[dr + 1] >= tc[0]:
              if time[dr] <= tc[1] and time[dr] >= tc[0]:
                 msd.append( ( 1e-3 * (space[dr + 1] - space[dr])**2 / (2 * (time[dr + 1] - time[dr]) ) ) )  #in um^2/ms
    ax23H.clear() 
    return msd    

def isd_calc(saddle_coordinates):
    if len(saddle_coordinates) > 1:
       isd = (saddle_coordinates[-1] - saddle_coordinates[0]) / (len(saddle_coordinates) - 1)
    else:
       isd = 0 
    return isd

def spot_boundary_spine(coordinate):
    proximity = False
    if abs(coordinate - start_spine) < 1e-6:
       proximity = True
    if abs(coordinate - end_spine) < 1e-6:
       proximity = True
    return proximity
        

def plot_corr(yratio, axis, axis_t, plot_flag):
  stim_patterns = ["stim"]
  parCount = 0
  spines = []
  spines_c = []
  spines_ic = []
  isd = []
  for stp in stim_patterns:
     tag_string = "NstimL_10.0_yratio_" + str(yratio) + "_spine_spacing_30.0_spacing_1.0_bfreq_20.0_coactive_1.0_mixed_True_trial_0" + stp + ".csv"
     parString = "NstimL_10.0_yratio_" + str(yratio) + "_spine_spacing_30.0_spacing_1.0_bfreq_20.0_coactive_1.0_mixed_True_trial_0"
     if stp == "stim":
        time_string = "timeC"
        space_string = "locC"
     else:    
        time_string = "timeIC"
        space_string = "locIC"
     avg_msd = []
     ic_isd_list = []
     c_isd_list = []
     print(tag_string)
     for tnum in range(len(trials)):
        if parCount < 1:
          phi = pd.read_csv("./" + parString + "phi.csv")
          sstart = pd.read_csv("./" + parString + "sstart.csv")
          final_cords = np.asarray(remove_nan( list(sstart.iloc[-1][1:]) ))
          print(tnum)
          print("Unfiltered: ", final_cords)
          ss = final_cords[np.where(final_cords > start_sync_range)]
          ss1 = ss[np.where(ss < end_sync_range)]
          c_isd_list.append( isd_calc(ss1) )
          spines_c.append(len(ss1))
          iss = list(final_cords[np.where(final_cords < start_sync_range)])
          print("Before removal: ", iss)
          for c in iss:
             proximal = spot_boundary_spine(c)
             if proximal:
               iss.remove(c)
          print("After: ", iss)
          iss_isd = isd_calc(iss)
          iss1 = list(final_cords[np.where(final_cords > end_sync_range)])
          print("Before removal: ", iss1)
          for c in iss1:
             proximal = spot_boundary_spine(c)
             if proximal:
               iss1.remove(c)
          print("After: ", iss1)
          iss1_isd = isd_calc(iss1)
          ic_isd = (iss_isd + iss1_isd) / 2.0
          ic_isd_list.append( ic_isd )
          spines_ic.append(len(iss) + len(iss1))
          finalN = len(remove_nan(list(phi.iloc[-1][1:])))
          curr_isd = (final_cords[-1] - final_cords[0]) / (finalN - 1)
          isd.append( curr_isd )
          spines.append(finalN)
        df = pd.read_csv("./" + tag_string)
        time = list( df[time_string] )
        space = np.asarray( df[space_string] ) * 1e6
        msd = MSD(space, time)
        avg_msd.append( sum(msd) / len(msd) )
        if tnum < len(trials) - 1:
           tag_string = tag_string.replace("trial_" + str(trials[tnum]), "trial_" + str(trials[tnum + 1]))
           parString = parString.replace("trial_" + str(trials[tnum]), "trial_" + str(trials[tnum + 1]))
     if plot_flag:
         axis_t.plot(trials, avg_msd, "b", marker = "*", linestyle = "--", label = "Stim Correlation")
     parCount +=1
#ax1H.plot(trials, spines, marker = "*")
     if plot_flag:
        axis.plot(trials, spines_c, "g", marker = "*", label = "Coherent zone N (10 $\mu m$)")
        axis.plot(trials, spines_ic, "m",  marker = "*", label = "Incoherant zone N (20 $\mu m$)")
        axis.plot(trials, np.asarray(spines_ic) + np.asarray(spines_c), "k",  marker = "*", label = "Incoherant zone N (20 $\mu m$)")
        axis.set_xlabel("Trials #")
        axis.set_ylabel("N #")
        axis_t.set_ylabel("Stim Correlation $\mu m^{2} / ms$")
        axis.legend(frameon = False, loc = "center left")
        axis_t.legend(frameon = False, loc = "upper left")
     return np.average(spines_c), np.average(spines_ic), spines_c, spines_ic, c_isd_list, ic_isd_list

avg_spines_c_list = []
avg_spines_ic_list = []
avg_total = []
yratio_list = [1.0]
axis_list = [ax11H, ax12H, ax21H]
axis_t_list = [ax11Ht, ax12Ht, ax21Ht]
for yratio in yratio_list:
   df = pd.DataFrame()
   axis_list[0].set_title("yBg/yStim = " + str(yratio))
   avg_spines_c, avg_spines_ic, spines_c, spines_ic, c_isd, ic_isd = plot_corr(yratio, axis_list[0], axis_t_list[0], True)
   axis_list.pop(0)
   axis_t_list.pop(0)
   avg_spines_c_list.append(avg_spines_c)
   avg_spines_ic_list.append(avg_spines_ic)
   avg_total.append(avg_spines_c + avg_spines_ic)
   df["spines_c"] = spines_c
   df["spines_ic"] = spines_ic
   df["IC_ISD"] = ic_isd
   df["C_ISD"] = c_isd
   df.to_csv(str(yratio) + "_" + "spine_stat.csv")
print(yratio_list)
print(avg_spines_c_list)
print(avg_spines_ic_list)
print(avg_total)
#ax21H.plot(yratio_list, avg_spines_c_list, label = "N Coherent")
#ax21H.plot(yratio_list, avg_spines_ic_list, label = "N InCoherent")
#ax21H.plot(yratio_list, avg_total, label = "Total")
#ax21H.set_xlabel("yBg/yStim")
#ax21H.legend(frameon = False)

tag_string = "NstimL_10.0_yratio_1.0_spine_spacing_30.0_spacing_1.0_bfreq_20.0_coactive_1.0_mixed_True_trial_0bg.csv"
avg_msd = []
for tnum in range(len(trials)):
   df = pd.read_csv("./" + tag_string)
   time_string = "timeIC"
   space_string = "locIC"
   ax1H.plot(df[time_string], [tnum] * len(df[time_string]), marker = "o", markerfacecolor = "b", markeredgecolor = "b", markersize = 1, linestyle = "") 
   if tnum < len(trials) - 1:
      tag_string = tag_string.replace("trial_" + str(trials[tnum]), "trial_" + str(trials[tnum + 1]))
ax1H.set_ylabel("Trial #")
ax1H.set_xlabel("Time (s)")

tag_string = "NstimL_10.0_yratio_1.0_spine_spacing_30.0_spacing_1.0_bfreq_20.0_coactive_1.0_mixed_True_trial_0stim.csv"
for tnum in range(len(trials)):
   df = pd.read_csv("./" + tag_string)
   time_string = "timeC"
   space_string = "locC"
   time = list( df[time_string] )
   space = np.asarray( df[space_string] ) * 1e6
   msd = MSD(space, time)
   avg_msd.append( sum(msd) / len(msd) )
   ax1H.plot(df[time_string], [21 + tnum] * len(df[time_string]), marker = "o", markerfacecolor = "g", markeredgecolor = "g", markersize = 1, linestyle = "") 
   if tnum < len(trials) - 1:
      tag_string = tag_string.replace("trial_" + str(trials[tnum]), "trial_" + str(trials[tnum + 1]))


plt.show() 

   

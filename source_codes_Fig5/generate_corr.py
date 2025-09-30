import tag
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("paper")

trials = np.arange(0, 15, 1)
#trials = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10]
#trials = [7, 8, 9, 10]
time_bin_size = 0.1

def remove_nan(lst):
    lst_nan = []
    for i in lst:
        if np.isnan(i) == False:
            lst_nan.append(i)
    return lst_nan        

def collated_db():
   tag_string = tag.tag_string
   tnum = 0
   collate_df = pd.DataFrame()
   for tr in trials:
      try: 
         print(tag_string)
         space_time = pd.read_csv("./" + tag_string + "bg.csv")
         coherent_times = space_time["timeC"]
         coherent_locs = space_time["locC"]

         ncoherent_times = space_time["timeIC"]
         ncoherent_locs = space_time["locIC"]


         all_locs = list(coherent_locs) + list(ncoherent_locs)
         all_times = list(coherent_times) + list(ncoherent_times)

         all_locs = np.round(np.asarray(all_locs) * 1e6, 1)
         df = pd.DataFrame()
         df["locs" + str(tnum)] = all_locs
         df["times" + str(tnum)] = all_times
         sorted_df = df.sort_values("locs" + str(tnum))
         sorted_df["indices" + str(tnum)] = np.arange(0, len(sorted_df),1)
         collate_df["locs" + str(tnum)] = sorted_df["locs" + str(tnum)]
         collate_df["times" + str(tnum)] = sorted_df["times" + str(tnum)]
         collate_df["indices" + str(tnum)] = sorted_df["indices" + str(tnum)]
         tag_string = tag_string.replace("trial_" + str(tnum), "trial_" + str(tnum + 1))
         tnum = tnum + 1
         collate_df.to_csv("./all_stim_data.csv")
      except:
         print("Trial: " + str(tnum) + " don't exist")
         tag_string = tag_string.replace("trial_" + str(tnum), "trial_" + str(tnum + 1))
         tnum = tnum + 1

def generate_cor(trial_num):
       df = pd.read_csv("./all_stim_data.csv")
       df_cor = pd.DataFrame()
       for d in [trial_num]:
             try:
                current_time_frame = list(df["times" + str(d)])
                current_locs = list(df["locs" + str(d)])
                tbins = np.arange(0, max(current_time_frame) + time_bin_size, time_bin_size)
                times_fired = []
                fired_locs = []
                proximities = []
                for tbin in range(len(tbins) - 1):
                     temp_locs = []
                     for current_time in range(len(current_time_frame)):
                        if current_time_frame[current_time] >= tbins[tbin]  and current_time_frame[current_time] <=tbins[tbin + 1]:
                            temp_locs.append(current_locs[current_time])
                     times_fired.append(len(temp_locs))
                     fired_locs.append(temp_locs)
                     if len(temp_locs) > 1:
                        proximities.append(1 / (max(temp_locs) - min(temp_locs)))
                     else:
                        proximities.append(np.inf)
                df_cor["fired_locs"] = fired_locs
                probs = np.asarray(times_fired) / sum(times_fired)
                df_cor["p"] = probs
                df_cor["proximity"] = proximities
                df_cor["fired_count"] = times_fired
             except:
                raise
                #print("trial: " + str(d) + "do not exist") 
       return df_cor 

def generate_entropy(pvals, distance_corr):
    sumE = 0
    scaled_sumE = 0
    print(pvals)
    for p in range(len(pvals)):
        if pvals[p] > 0:
           sumE = sumE + -1 * pvals[p] * np.log(pvals[p])
           scaled_sumE = scaled_sumE + -1 * distance_corr[p] * pvals[p] * np.log(pvals[p])
    print(sumE)
    return sumE, scaled_sumE
        

collated_db()

trials = np.arange(0, 15, 1)
df_trials = pd.DataFrame()
n_list = []
e_list = []
scaled_e_list = []
isd_list = []
tag_string = tag.tag_string
for tnum in trials:
   try:
     df_e = generate_cor(tnum)
     entropy, scaled_entropy = generate_entropy(df_e["p"], df_e["proximity"])
     phiFile = pd.read_csv("./" + tag_string + "phi.csv")
     cordFile = pd.read_csv("./" + tag_string + "sstart.csv")
     if len(phiFile) > 0:
        num_spines = 0
        for n in list(phiFile.iloc[-1][1:]):
           if np.isnan(n):
             num_spines = num_spines + 0
           else:
             num_spines = num_spines + 1
        saddle_starts = remove_nan(cordFile.iloc[-1][1:])
        isd = (saddle_starts[-1] - saddle_starts[0]) / (num_spines - 1)
     else:
        num_spines = 0
        isd = np.nan
     n_list.append(num_spines)
     e_list.append(entropy)
     scaled_e_list.append(scaled_entropy)
     isd_list.append(isd)
     print(entropy)
   except:
     #raise
     print("Trial: " + str(tnum) + " don't exist")
   tag_string = tag_string.replace("trial_" + str(tnum), "trial_" + str(tnum + 1))

df_trials["N"] = n_list
df_trials["Entropy"] = e_list
df_trials["scaled Entropy"] = scaled_e_list
df_trials["ISD"] = np.asarray(isd_list) * 1e6

df_trials.to_csv("./entropy_across_trials.csv")
print(df_trials)
df = pd.read_csv("./entropy_across_trials.csv")
sorted_en = df.sort_values("Entropy")
plt.plot(sorted_en["Entropy"], sorted_en["N"], "o")
plt.ylabel("N")
plt.xlabel("Entropy")
plt.ylim([0,5])
plt.title("F: " + str(tag.bfreq))
#plt.plot(sorted_en["Entropy"], sorted_en["ISD"], marker = "*", linestyle = None)
plt.show()

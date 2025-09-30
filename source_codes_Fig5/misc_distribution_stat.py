import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("poster")




xlimit = 80

time_bins = np.arange(0, 83.0, 1.0)

def make_tag_S(bfreq, kcat):
    file_q = "num_spines_S_" + str(kcat) + "_" + str(bfreq) + ".csv"
    file_t = "num_spines_time_S_" + str(kcat) + "_" + str(bfreq) + ".csv"
    return file_q, file_t

def make_tag(bfreq, kcat):
    file_q = "num_spines_" + str(kcat) + "_" + str(bfreq) + ".csv"
    file_t = "num_spines_time_" + str(kcat) + "_" + str(bfreq) + ".csv"
    return file_q, file_t

def make_tag_AS(bfreq, kcat):
    file_q = "num_spines_AS_" + str(kcat) + "_" + str(bfreq) + ".csv"
    file_t = "num_spines_time_AS_" + str(kcat) + "_" + str(bfreq) + ".csv"
    return file_q, file_t

def make_trial_list(file_df):
    trials = list(file_df.columns[1:])  
    return trials

def count_bins(df_time, df_num):
    bin_count = [0] * len(time_bins)
    for i in range(len(df_time)):
       for j in range(len(time_bins) - 1):
          if df_time.iloc[i] >= time_bins[j] and df_time.iloc[i] <= time_bins[j + 1]:
             bin_count[j] = df_num.iloc[i]
             break
    return bin_count

def plot_trials(df_q, df_t, axis):
    trial_list = make_trial_list(df_q)
    for trial in trial_list:
        counts = count_bins(df_t[str(trial)], df_q[str(trial)])
        axis.plot(time_bins, counts, label = trial)
    axis.set_xlim(0, 80)
    axis.legend()

def plot_avg(df_q, df_t, axis):
    trial_list = make_trial_list(df_q)
    count_sum = np.asarray([0] * len(time_bins))
    for trial in trial_list:
        count_sum = count_sum + np.asarray(count_bins(df_t[str(trial)], df_q[str(trial)]))
    avg = np.asarray(count_sum) / len(trial_list)
    axis.plot(time_bins, avg, label = "Avg")
    axis.set_xlim(0, 80)
    axis.legend()

def plot_avg_density(df_q, df_t, axis, length, label):
    trial_list = make_trial_list(df_q)
    count_sum = np.asarray([0] * len(time_bins))
    for trial in trial_list:
        count_sum = count_sum + np.asarray(count_bins(df_t[str(trial)], df_q[str(trial)]))
    avg_count = np.asarray(count_sum) / len(trial_list)
    avg_density = avg_count / length
    axis.plot(time_bins, avg_density, label = "Avg. " + label)
    axis.set_xlim(0, 80)
    axis.legend()
    return time_bins, avg_density

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

kcat_list = [20.0, 50.0, 200.0]
axes = [ax1, ax2, ax3]

for kcat in kcat_list:
   Sbfreq = 10.0
   file_q, file_t = make_tag_S(Sbfreq, kcat)
   df_num = pd.read_csv(file_q)
   df_t = pd.read_csv(file_t)
   time,density = plot_avg_density(df_num, df_t, axes[0], 10, "S")
   print(df_num)
   df = pd.DataFrame()
   df["time"] = time
   df["density"] = density
   df.to_csv("./time_density_" + str(kcat) + "_S" + str(Sbfreq) + ".csv")

   ASbfreq = 2 * Sbfreq
   file_q, file_t = make_tag_AS(ASbfreq, kcat)
   df_num = pd.read_csv(file_q)
   df_t = pd.read_csv(file_t)
   time, density = plot_avg_density(df_num, df_t, axes[0], 20, "AS")
   axes[0].set_title("kcat: " + str(kcat) + ", numSamples: " + str(len(df_num.columns[1:])))
   df = pd.DataFrame()
   df["time"] = time
   df["density"] = density
   df.to_csv("./time_density_" + str(kcat) + "_AS" + str(ASbfreq) + ".csv")
   if axes[0] == ax2:
      axes[0].set_ylabel("Spine density #/$\mu m$")
   axes.pop(0)

plt.show()


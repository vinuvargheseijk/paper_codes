import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
import xml.etree.ElementTree as ET
import readXML
import seaborn as sns
sns.set_context("poster")

dt = 0.2
runtime = 85
max_index = int(runtime / 0.2)
start_sync = 24e-6
end_sync = 36e-6
sync_length = 10 #\mu m
async_length = 20
boundary_spine1 = 11e-6
boundary_spine2 = 39e-6
stim_time = 80
Length = 50e-6
dendR = 0.5e-6
diffL = 20e-9
v_cyt = np.pi * dendR**2 * Length
Na = 6.0223e23


def remove_nan(list_v):
    nonNan = []
    for i in list_v:
         if np.isnan(i) == False:
            nonNan.append(i)
    return nonNan

def remove_boundary_spines(list_v):
    nonB = []
    for i in list_v:
       if i > boundary_spine1 and i < boundary_spine2:
           nonB.append(i)
    return nonB

def diff_2_points(data):
    diff_l = []
    for i in range(len(data) - 1):
        diff_l.append(abs(data[i + 1] - data[i]))
    return diff_l


def synchronous_count(directory, tag_string):
    coord_file = "sstart.csv"
    df_coord = pd.read_csv(directory + "/" + tag_string + "sstart.csv")
    df_time = pd.read_csv(directory + "/" + tag_string + "shapeTimes.csv")
    num_spines_in_time = []
    time = []
    for l in range(len(df_coord)):
       nonNan_list = remove_nan(list(df_coord.iloc[l][1:]))
       count = 0
       for ns in range(len(nonNan_list)):
           if nonNan_list[ns] > start_sync and nonNan_list[ns] < end_sync:
              count += 1
       num_spines_in_time.append(count)
       time.append(df_time.iloc[l]["0"])
    max_saved_length = len(num_spines_in_time)
    if max_saved_length < max_index:
       fill_num = []
       fill_time = []
       for extra in range(max_index - max_saved_length):
           fill_num.append(np.nan)
           fill_time.append(time[-1] + dt * extra)
       num_spines_in_time = num_spines_in_time + fill_num
       time = time + fill_time
    return time, num_spines_in_time   

def asynchronous_count(directory, tag_string):
    coord_file = "sstart.csv"
    df_coord = pd.read_csv(directory + "/" + tag_string + "sstart.csv")
    df_time = pd.read_csv(directory + "/" + tag_string + "shapeTimes.csv")
    num_spines_in_time = []
    time = []
    for l in range(len(df_coord)):
       nonNan_list = remove_nan(list(df_coord.iloc[l][1:]))
       count = 0
       for ns in range(len(nonNan_list)):
           if nonNan_list[ns] < start_sync or nonNan_list[ns] > end_sync:
              if nonNan_list[ns] > boundary_spine1 and nonNan_list[ns] < boundary_spine2:
                 count += 1
       num_spines_in_time.append(count)
       time.append(df_time.iloc[l]["0"])
    max_saved_length = len(num_spines_in_time)
    if max_saved_length < max_index:
       fill_num = []
       fill_time = []
       for extra in range(max_index - max_saved_length):
           fill_num.append(np.nan)
           fill_time.append(time[-1] + 0.2 * extra)
       num_spines_in_time = num_spines_in_time + fill_num
       time = time + fill_time
    return time, num_spines_in_time   

def count_spines(directory, tag_string):
    phi_file = "phi.csv"
    df_phi = pd.read_csv(directory + "/" + tag_string + "phi.csv")
    df_time = pd.read_csv(directory + "/" + tag_string + "shapeTimes.csv")
    df_coord = pd.read_csv(directory + "/" + tag_string + "sstart.csv")
    num_spines_in_time = []
    time = []
    for l in range(len(df_phi)):
        nonNan_coords = remove_nan(list(df_coord.iloc[l][1:]))
        nonB_spines = remove_boundary_spines(nonNan_coords)
        num_spines_in_time.append(len(nonB_spines))
        time.append(df_time.iloc[l]["0"])
    max_saved_length = len(num_spines_in_time)
    if max_saved_length < max_index:
       fill_num = []
       fill_time = []
       for extra in range(max_index - max_saved_length):
           fill_num.append(np.nan)
           fill_time.append(time[-1] + dt * extra)
       num_spines_in_time = num_spines_in_time + fill_num
       time = time + fill_time
    return time, num_spines_in_time   

def find_max_spines_across_trials(df):
    max_spines = 0
    for i in range(len(df)):
       for trial in found_trials:
         if df.iloc[i][str(trial)] > max_spines:
            max_spines = df.iloc[i][str(trial)]
    return max_spines

def calc_cytConc(directory, tag_string):
     filename = tag_string+"cytosol.xml"
     tree = ET.parse(directory + "/" + filename)
     sum_an,max_an_y,an_y=readXML.plotXML(filename,[tree])
     num_voxels = Length / diffL
     avg_conc = 1e3 * np.asarray(sum_an) / (v_cyt * Na)
     chemTimes = pd.read_csv(tag_string + "chemTimes.csv")
     chemTimes = list(chemTimes["time"])
     print("Lengths: ", len(chemTimes), len(avg_conc))
     return chemTimes, avg_conc
     
def calc_Nmolecules(directory, tag_string, molecule):
     filename = tag_string + str(molecule) + ".xml"
     tree = ET.parse(directory + "/" + filename)
     sum_an,max_an_y,an_y=readXML.plotXML(filename,[tree])
     num_voxels = Length / diffL
     #avg_conc = 1e3 * np.asarray(sum_an) / (v_cyt * Na)
     total_number = np.asarray(sum_an)
     chemTimes = pd.read_csv(tag_string + "chemTimes.csv")
     chemTimes = list(chemTimes["time"])
     print("Lengths: ", len(chemTimes), len(avg_conc))
     return chemTimes, total_number
         
def avg_count(filenames):
    time = pd.read_csv(filenames[0])
    spine_count = pd.read_csv(filenames[1])
    max_spines = find_max_spines_across_trials(spine_count)
    avg_spines = []
    for ns in range(len(spine_count)):
        nonNans = remove_nan(list(spine_count.iloc[ns][1:]))
        avg_spines.append(sum(nonNans) / len(found_trials))
    return time[str(found_trials[0])], avg_spines

def avg_cytConc(filenames):
    time = pd.read_csv(filenames[0])
    cytConc = pd.read_csv(fielnames[1])
    

def add_hatch(axis, scale):
   stim_hatch_axis = Rectangle((0, 0), stim_time, scale * 1)
   axis.add_patch(stim_hatch_axis)
   axis.set_xlim(0, stim_time)


def make_tag(trial, kcat):
    tag_string = "NstimL_10.0_kcat_enz_tiam_" + str(kcat) + "_spine_spacing_30.0_spacing_1.0_bfreq_" + str(bfreq) + "_coactive_1.0_mixed_True_trial_" + str(trial)    
    return tag_string

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
trial_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
kcat_list = [20.0, 50.0, 200.0]


bfreq = 20.0
directory = "./"
label_count = 0
df_S_density = pd.DataFrame()
df_S_density_time = pd.DataFrame()
df_AS_density = pd.DataFrame()
df_AS_density_time = pd.DataFrame()
for kcat in kcat_list:
   found_trials = []
   df_count = pd.DataFrame()
   df_time = pd.DataFrame()
   df_count_S = pd.DataFrame()
   df_time_S = pd.DataFrame()
   df_count_AS = pd.DataFrame()
   df_time_AS = pd.DataFrame()
   for trial in trial_list:
       tag_string = make_tag(trial, kcat)
       print(tag_string)
       try:
         time, num_spine_list = count_spines(directory, tag_string)
         time_S, num_spine_S = synchronous_count(directory, tag_string)
         time_AS, num_spine_AS = asynchronous_count(directory, tag_string)
         chemTimes, avg_conc = calc_cytConc(directory, tag_string)
         ax4.plot(chemTimes, avg_conc, label = "cytosol Conc.")
         chemTimes, number_of_molecules = calc_Nmolecules(directory, tag_string, "cytosol")
         ax5.plot(chemTimes, number_of_molecules, label = "Cytosol")
         chemTimes, number_of_molecules = calc_Nmolecules(directory, tag_string, "membrane")
         ax6.plot(1e3 * np.array(diff_2_points(number_of_molecules)) / (v_cyt * Na) )
         ax6.set_title("Dumped cytosol conc")
         ax5.plot(chemTimes, number_of_molecules, label = "Membrane")
         if label_count == 0:
            ax4.plot(chemTimes, [0.4] * len(avg_conc), 'k', label = "Threshold")
         df_count_S[str(trial)] = num_spine_S
         df_time_S[str(trial)] = time_S
         df_count_AS[str(trial)] = num_spine_AS
         df_time_AS[str(trial)] = time_AS
         df_count[str(trial)] = num_spine_list
         df_time[str(trial)] = time
         if label_count == 0:
            label = "Trials"
         else:
            label = None
         ax1.plot(time_AS, num_spine_AS, "grey", label = label)
         ax2.plot(time_S, num_spine_S, "grey", label = label)
         label_count += 1    
         found_trials.append(trial)
       except Exception as exception:
         print("Exception: ", exception)
   print(df_time)
   df_count.to_csv("./num_spines_" + str(kcat) + "_" + str(bfreq) + ".csv")
   df_time.to_csv("./num_spines_time_" + str(kcat) + "_" + str(bfreq) +".csv")
   df_count_S.to_csv("./num_spines_S_" + str(kcat) + "_" + str(bfreq) + ".csv")
   df_time_S.to_csv("./num_spines_time_S_" + str(kcat) + "_" + str(bfreq) + ".csv")
   df_count_AS.to_csv("./num_spines_AS_" + str(kcat) + "_" + str(bfreq) +".csv")
   df_time_AS.to_csv("./num_spines_time_AS_" + str(kcat) + "_" + str(bfreq) + ".csv")

plt.tight_layout()
plt.show()

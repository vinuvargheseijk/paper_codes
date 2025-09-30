from functools import reduce
import sys
import argparse
import numpy as np
import pandas as pd



parser = argparse.ArgumentParser()
parser.add_argument('--kagg', type = float)
parser.add_argument('--ns', type = float)
parser.add_argument('--mu0', type = float)
parser.add_argument('--dt', type = float)
parser.add_argument('--spacing', type = float)
parser.add_argument('--cytD', type = float)
parser.add_argument('--phi_entire', type = float)
############################################################
parser.add_argument('--ypos', type = str, required  = True)
parser.add_argument('--freq', type = float, required  = True)
parser.add_argument('--freq_ratio', type = float, required  = True)
parser.add_argument('--bfreq', type = str, required  = True)
parser.add_argument('--bs', type = str, required  = True)
parser.add_argument('--nb', type = str, required  = True)
parser.add_argument('--CaTau', type = float, required  = True)
parser.add_argument('--s', type = float, required  = True)
#############################################################
args = parser.parse_args()
Length = 10
Kb = 40.0
runtime = 3
bigDia = 1.0
smallDia = 1.0
comptLen = 0.01
CaScale = 1.0
#peclet = args.pe


def calcTaper():
  taper = (bigDia - smallDia) / (Length / comptLen)
  return taper
phi_entire = args.phi_entire
Rac_Metab_Reac_Kf = 0.0
Tiam_Metab_Reac_Kf = 0.0

Enz_Tiam_kcat = 50.0
#Enz_tot_CaM_CaMKII_kcat = 4.0

Rac_Rev_Reac_Kf = 0.5
Tiam_Rev_Reac_Kf = 1.0

Enz_tot_CaM_CaMKII_Km = 0.05
Enz_Tiam_Km = 0.08
CaTau = args.CaTau
mem_diffConst = 0.0
cyt_diffConst =  args.cytD * 1e-12
freq = args.freq
dtb = 1/args.freq
freq_ratio = args.freq_ratio
freq2 = args.freq_ratio * args.freq
dtb_2nd_spine = 1/ freq2
dtb_list = [dtb, dtb_2nd_spine]
num_bursts = int(args.nb)
end_stim_time = num_bursts * dtb + 2.0
firing_rate = freq
############################################################3
bg = False
bg_int = float(args.bs)
bg_freq = float(args.bfreq)
#bg_start = end_stim_time
bg_start = 0.0
print("BG start: ", bg_start)
bg_end = runtime
bg_dt = 1 / bg_freq
#num_bg_times = (bg_end - bg_start) / bg_dt
num_bg_times = num_bursts
np.random.seed(3)
bg_locs = np.asarray( np.random.random(int(num_bg_times)) ) * 10 * 1e-6
print("BG locaions: ", bg_locs)
y_pos = float(args.ypos)
Source = float(args.s)
###########################################################
#cyt_motor_rate = cyt_diffConst * peclet / (bigDia * 1e-6)
#print("motor rate: ", cyt_motor_rate)
num_spines = int(args.ns)   #num_spines are num_clusters here
cluster_size = 1
cluster_spacing = 1.0 * 1e-6
scalF = 1.0
mu0_exp = args.mu0
k_agg =  args.kagg
k_agg2 = -k_agg / 2.0
spacing = args.spacing * 1e-6
fix_conc = False
start_spine = ( Length * 1e-6 - (num_spines - 1) * spacing ) / 2.0  #Check if this is not zero, which means the length is not enought to support num_spines.
#start_spine = 2 * 1e-6  #Check if this is not zero, which means the length is not enought to support num_spines.
print("START SPINE: ", start_spine)
p_code = 2
taper = calcTaper()
stim_locations = []
stim_times = []
###########################################################
glut_stim_times = []
glut_stim_locations = []
glut_stim_y = []
delay = 0
stim_type = []
burst_num_pulse = 10
ypos_ratio = 1.0
#########
########
try:
    dfTimes = pd.read_csv("./expTimeF" + str(freq) + ".csv")
    exp_times = list(dfTimes["times"])
except:
    print("File not here")
    dfTimes = pd.DataFrame()
    exp_times = np.random.exponential(scale = 1 / firing_rate, size = num_spines * num_bursts)
    dfTimes["times"] = exp_times
    dfTimes.to_csv("./expTimeF" + str(freq) + ".csv")
for i in range(num_spines * num_bursts):
    stim_type.append({"stim" : "s"})
bg_stim_type = "s"
#firing_time = exp_times[0]
firing_time = 0.5

for i in range(num_bursts):
   for n in range(num_spines):
       glut_stim_times.append( 0.5 + i * dtb )
       #glut_stim_times.append( delay + i * dtb )
       glut_stim_locations.append(start_spine + n * spacing)
       position = y_pos - n * ypos_ratio * y_pos
       if position > 0:
          glut_stim_y.append(position)
       else:   
          glut_stim_y.append(y_pos)
       #########
       #firing_time = firing_time + exp_times[i * n]
       #########
##########################################################   
print("Freq and poisson rate: ", freq, num_bursts / sum(exp_times))
#glut_stim_locations.append(glut_stim_locations[-1] + 0.2e-6)
#glut_stim_times.append(glut_stim_times[-1])
#glut_stim_y.append(y_pos)
print("STIM LOCATIONS: ", glut_stim_locations)
print("STIM TIMES: ", glut_stim_times)
print("STIM Y: ", glut_stim_y)
if len(stim_locations) == 1:
    listToStr = 'One'
tag_string = "Source_" + str(args.s) + "_nb_" + str(num_bursts) + "_ns_" + str(num_spines) + "_ypos_" + str(y_pos) + "_spacing_" + str(spacing) + "_freq_" + str(args.freq) + "_CaTau_" + str(CaTau) + "_phi_entire_" + str(phi_entire)
print("Stim locations: ", stim_locations)

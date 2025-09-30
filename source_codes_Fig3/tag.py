import argparse
import numpy as np


chemTurnOff = True
Length = 20
Kb = 40.0
bigDia = 1.0
smallDia = 1.0
comptLen = 0.1
CaScale = 1.0
#peclet = args.pe
init_dummy = True

parser = argparse.ArgumentParser()
parser.add_argument('--kagg', type = float)
parser.add_argument('--ns', type = float)
parser.add_argument('--mu0', type = float)
parser.add_argument('--dt', type = float)
parser.add_argument('--spacing', type = float)
parser.add_argument('--cytD', type = float)
parser.add_argument('--cyt_conc', type = float)
parser.add_argument('--dummy_spine_spacing', type = float, required  = True)
parser.add_argument('--spine_create_delay', type = float, required  = True)
############################################################
if chemTurnOff == False:
   parser.add_argument('--typeStim', type = str)
   parser.add_argument('--yposR', type = float, required  = True)
   parser.add_argument('--ypos', type = float, required  = True)
   parser.add_argument('--freq', type = float, required  = True)
   parser.add_argument('--bfreq', type = str, required  = True)
   parser.add_argument('--bs', type = str, required  = True)
   parser.add_argument('--nb', type = str, required  = True)
   parser.add_argument('--s', type = float, required  = True)
   parser.add_argument('--NstimL', type = float, required  = True)
   parser.add_argument('--TK', type = int, required  = True)
   parser.add_argument('--trial', type = int, required  = True)
args = parser.parse_args()

if chemTurnOff:
   cyt_conc = args.cyt_conc
else:
   cyt_conc = 0.00005  #Basal conc of activated IRSp53

CaTau = 0.08
if chemTurnOff == False:
   spacing = args.spacing * 1e-6
   Rac_Metab_Reac_Kf = 0.0
   Tiam_Metab_Reac_Kf = 0.0
   TK = args.TK
   chem_perturb = False
   Enz_Tiam_kcat = 50.0
   #Enz_tot_CaM_CaMKII_kcat = 4.0

   Rac_Rev_Reac_Kf = 0.5
   Tiam_Rev_Reac_Kf = 1.0

   Enz_tot_CaM_CaMKII_Km = 0.05
   Enz_Tiam_Km = 0.08
   freq = args.freq
   dtb = 1/args.freq
   num_bursts = int(args.nb)
   end_stim_time = num_bursts * dtb + 2.0
   firing_rate = freq
   Source = float(args.s)
   tiam_dependency = False
   bg = False
   num_stim_locations = args.NstimL
   start_stim = ( Length * 1e-6 - (num_stim_locations - 1) * spacing ) / 2.0  #Check if this is not zero, which means the length is not enought to support num_spines.
   glut_stim_times = []
   glut_stim_locations = []
   glut_stim_y = []
   delay = 0
   stim_type = []
   burst_num_pulse = 10
   ypos_ratio = 1.0
   bg_stim_type = args.typeStim


settling_time = 1.0
spine_create_delay = args.spine_create_delay
def calcTaper():
  taper = (bigDia - smallDia) / (Length / comptLen)
  return taper


dimer_conc = 0.0025

trial_plots = True

mem_diffConst = 0.0
cyt_diffConst =  args.cytD * 1e-12
dummy_spine_spacing = args.dummy_spine_spacing * 1e-6
###########################################################
num_spines = int(args.ns)   #num_spines are num_clusters here
mu0_exp = args.mu0
k_agg =  args.kagg
k_agg2 = -k_agg / 2.0
fix_conc = False
start_spine = ( Length * 1e-6 - (num_spines - 1) * dummy_spine_spacing ) / 2.0  #Check if this is not zero, which means the length is not enought to support num_spines.
taper = calcTaper()
###########################################################
def stim_loc_specific(stim_loc, num_bursts, dtb_loc, offset, release_height):
   for i in range(int(num_bursts)):
       glut_stim_times.append( settling_time + i * dtb_loc + offset )  #The stimuli for second spine is offset till the stimuli to the first is over.
       glut_stim_locations.append( stim_loc )
       glut_stim_y.append( release_height )

if chemTurnOff:
   init_mature_spine = True
else:
   init_mature_spine = False
init_mature_spine_index = 0
dummy_spine_locations = []
spine_create_times = []
for n in range(0, num_spines):
    dummy_spine_locations.append(start_spine + n * dummy_spine_spacing)
    spine_create_times.append(settling_time + n * args.spine_create_delay)

deletion_delay = 0
if chemTurnOff == False:
   stim_loc_list = []
   for nstim in range(int(num_stim_locations)):
       stim_loc_list.append(start_stim + nstim * spacing)

   stim_loc_specific(stim_loc_list[0], int(num_bursts), dtb, args.spine_create_delay, args.ypos * 1e-6)
   #stim_loc_specific(stim_loc_list[2], int(num_bursts / 4), dtb, args.spine_create_delay, args.ypos * 1e-6)
   #stim_loc_specific(stim_loc_list[3], int(num_bursts / 4), dtb, args.spine_create_delay, args.ypos * 1e-6)
   #stim_loc_specific(stim_loc_list[4], int(num_bursts / 4), dtb, args.spine_create_delay, args.ypos * 1e-6)

   runtime = glut_stim_times[-1] + deletion_delay + 20
   print("GLUT STIM LOCATIONS: ", glut_stim_locations)
   print("GLUT STIM times: ", glut_stim_times)
   stim_start_time = min(glut_stim_times) - 1.0
   stim_end_time = max(glut_stim_times) + 1.0
   print("Stim start and end times: ", stim_start_time, stim_end_time)
   print("STIM TIMES: ", glut_stim_times)
   print("STIM Y: ", glut_stim_y)
   for i in range(int(len(glut_stim_times))):
        stim_type.append({"stim" : args.typeStim})
else:
   runtime = args.spine_create_delay + 60

print("Spine locations: ", dummy_spine_locations)
print("Spine create times: ", spine_create_times)
##########################################################   





if chemTurnOff:
     tag_string = "CytConc_" + str(args.cyt_conc) + "_ns_" + str(num_spines) + "_spine_create_delay_" + str(args.spine_create_delay) + "_spine_spacing_" + str(args.dummy_spine_spacing)
else:
     tag_string = "NstimL_" + str(args.NstimL) + "_nb_" + str(num_bursts) + "_ns_" + str(num_spines) + "_spacing_" + str(args.spacing) + "_freq_" + str(args.freq) + "_ypos_" + str(args.ypos) + "_spine_create_delay_" + str(args.spine_create_delay) + "_spine_spacing_" + str(args.dummy_spine_spacing) + "_TK_" + str(args.TK) + "_trial_" + str(args.trial)

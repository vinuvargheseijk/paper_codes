# Source Generated with Decompyle++
# File: tag.cpython-310.pyc (Python 3.10)

import argparse
import numpy as np
import pandas as pd
import random
import poisson_time as PT
chemTurnOff = False
Length = 50
Kb = 40
bigDia = 1
smallDia = 1
comptLen = 0.1
CaScale = 1
init_dummy = True
boundary_padding = 10e-06
parser = argparse.ArgumentParser()
parser.add_argument('--kagg', type = float)
parser.add_argument('--ns', type = float)
parser.add_argument('--mu0', type = float)
parser.add_argument('--dt', type = float)
parser.add_argument('--spacing', type = float)
parser.add_argument('--cyt_conc', type = float)
parser.add_argument('--spine_create_delay', type = float, required = True)
if chemTurnOff == False:
    parser.add_argument('--typeStim', type = str)
    parser.add_argument('--yposR', type = float, required = True)
    parser.add_argument('--ypos', type = float, required = True)
    parser.add_argument('--freq', type = float, required = True)
    parser.add_argument('--bfreq', type = float, required = True)
    parser.add_argument('--kcat_enz_tiam', type = float, required = True)
    parser.add_argument('--nb', type = str, required = True)
    parser.add_argument('--s', type = float, required = True)
    parser.add_argument('--NstimL', type = float, required = True)
    parser.add_argument('--dummy_spine_spacing', type = float, required = True)
    parser.add_argument('--coactive', type = float, required = True)
    parser.add_argument('--trial', type = int, required = True)
    parser.add_argument('--bgy', type = int, required = True)
    parser.add_argument('--mixed', type = bool, required = True)
runtime = 85
args = parser.parse_args()
ypos = args.ypos
num_spines = args.ns
dummy_spine_locations = []
spine_create_times = []
cyt_conc = args.cyt_conc
start_spine = (Length * 1e-06 - (num_spines - 1) * args.dummy_spine_spacing * 1e-06) / 2
settling_time = 1
num_spines = int(args.ns)
for n in range(0, num_spines):
    dummy_spine_locations.append(start_spine + n * args.dummy_spine_spacing * 1e-06)
    spine_create_times.append(settling_time + args.spine_create_delay)
print('SPINE LOCATIONS: ', dummy_spine_locations)
print('SPINE CREATION TIME: ', spine_create_times)
spine_create_delay = args.spine_create_delay
CaTau = 0.08
if chemTurnOff == False:
    dt_b = 0.025
    kcat_enz_tiam = args.kcat_enz_tiam
    spacing = args.spacing * 1e-06
    Rac_Metab_Reac_Kf = 0
    Tiam_Metab_Reac_Kf = 0
    chem_perturb = False
    Enz_Tiam_kcat = 50
    Rac_Rev_Reac_Kf = 0.5
    Tiam_Rev_Reac_Kf = 1
    Enz_tot_CaM_CaMKII_Km = 0.05
    Enz_Tiam_Km = 0.08
    freq = args.freq
    dtb = 1 / args.freq
    num_bursts = int(args.nb)
    end_stim_time = num_bursts * dtb + 2
    firing_rate = freq
    Source = float(args.s)
    tiam_dependency = True
    bg = True
    num_stim_locations = args.NstimL
    if bg:
        bg_spacing = spacing
        #start_bg_stim = (Length * 1e-06 - (num_bg_locations - 1) * bg_spacing) / 2
        start_bg_stim = start_spine
        bg_stim_type = "b"
        bg_ypos = args.bgy * ypos * 1e-06
        bg_stim_y_ratio = bg_ypos / (ypos * 1e-6)
    #start_stim = (Length * 1e-06 - (num_stim_locations - 1) * spacing) / 2
    #start_stim = dummy_spine_locations[-1] - (num_stim_locations - 1) * spacing
    start_stim = 25e-6
    glut_stim_times = []
    glut_stim_locations = []
    glut_stim_y = []
    delay = 0
    stim_type = []
    burst_num_pulse = 10
    ypos_ratio = 1
settling_time = 0
spine_create_delay = args.spine_create_delay

def calcTaper():
    taper = (bigDia - smallDia) / Length / comptLen
    return taper

dimer_conc = 0.0025 - cyt_conc
mem_diffConst = 0
mu0_exp = args.mu0
k_agg = args.kagg
k_agg2 = -k_agg / 2
fix_conc = False
taper = calcTaper()

def stim_loc_specific(stim_loc, num_bursts, dtb_loc, offset, release_height):
    for i in range(int(num_bursts)):
        print('ADDING STIM: ', stim_loc)
        glut_stim_times.append(settling_time + i * dtb_loc + offset)
        glut_stim_locations.append(stim_loc)
        glut_stim_y.append(release_height)

df_stim_times = pd.DataFrame()

def poisson_stim_loc_specific(stim_loc, num_bursts, offset, release_height):
    exp_times = np.random.exponential(1 / freq, num_bursts)
    previous = 0
    individual_stim_list = []
    for i in range(int(num_bursts)):
        glut_stim_times.append(settling_time + previous + exp_times[i] + offset)
        individual_stim_list.append(settling_time + previous + exp_times[i] + offset)
        glut_stim_locations.append(stim_loc)
        glut_stim_y.append(release_height)
        previous = previous + exp_times[i]
    df_stim_times[str(stim_loc * 1e+06)] = individual_stim_list
    print(glut_stim_times)

bfreq = args.bfreq

def poisson_bg(bg_sources, coherence):
    print("BG sources: ", bg_sources)
    duration = runtime - 5
    dt = 0.001
    #scaling_ratio = PT.scaling_bg_freq(args.bfreq, duration, dt)
    scaling_ratio = 3.0
    print("Scaling ratio: ", scaling_ratio)
    if args.mixed:
       print("Mixed stimulii")
       if coherence:
          times = PT.IpSpike(args.bfreq, duration, dt)
       else:
          times = PT.pSpike( ( 1 / scaling_ratio ) * args.bfreq, duration)
    else:
       print("Homogenous stimuli everywhere")
       times = PT.pSpike(args.bfreq, duration)
    bg_times = []
    rand_locs = []
    for et in times:
        temp_loc_index = random.randint(0, len(bg_sources) - 1)
        rand_locs.append(bg_sources[temp_loc_index])
        bg_times.append(et)
    return (bg_times, rand_locs)    

allowDeletion = True
deletion_delay = 0
if chemTurnOff == False:
    stim_loc_list = []
    if int(args.NstimL) > 0:
        for nstim in range(int(num_stim_locations)):
            stim_loc_list.append(start_stim + nstim * spacing)
        print('STIM LOCS: ', stim_loc_list)
    bsflag1 = []
    blflag1 = []
    byflag1 = []

    def arrange_stim(bg_locs, bg_times, bg_ypos):
        for bg_ti in range(len(bg_times)):
            if bg_stim_type == 's':
                bsflag1.append(bg_times[bg_ti])
                blflag1.append(bg_locs[bg_ti])
                byflag1.append(bg_ypos)
            if bg_stim_type == 'b':
                single_burst_time = []
                single_burst_loc = []
                single_burst_y = []
                for i in range(burst_num_pulse):
                    single_burst_time.append(bg_times[bg_ti] + i * dt_b)
                    single_burst_loc.append(bg_locs[bg_ti])
                    single_burst_y.append(bg_ypos)
                bsflag1.append(single_burst_time)
                blflag1.append(single_burst_loc)
                byflag1.append(single_burst_y)
 
    if bg:
        df_bg_C = pd.DataFrame()
        df_bg_IC = pd.DataFrame()
        
        bg_sources = []
        bg_start = settling_time
        bg_pos = 0
        bg_l = 0
        while bg_pos < (start_stim - bg_spacing):
            print(bg_pos)
            bg_sources.append(start_bg_stim + bg_l * bg_spacing)
            bg_pos = bg_sources[-1]
            bg_l = bg_l + 1
        print("Stims: ", stim_loc_list)
        bg_pos = 0
        bg_l = 1
        while bg_pos < dummy_spine_locations[-1]:
            bg_sources.append(stim_loc_list[-1] + bg_l * bg_spacing)
            bg_pos = bg_sources[-1]
            bg_l = bg_l + 1
        print("Sources: ", bg_sources)  
        print("Len of Sources: ", len(bg_sources))  
 
        (bg_times, bg_locs) = poisson_bg( bg_sources, False )
        print("Time of release: ", bg_times)
        (stim_times, stim_locs) = poisson_bg( stim_loc_list, True )
   
        diff_count = len(bg_times) - len(stim_times)
        if diff_count > 0:
            for dc in range(abs(diff_count)):
                bg_times.pop(-1)
                bg_locs.pop(-1)
        if diff_count < 0:
            for dc in range(abs(diff_count)):
                stim_times.pop(-1)
                stim_locs.pop(-1)
        arrange_stim(stim_locs, stim_times, args.ypos * 1e-6)
        arrange_stim(bg_locs, bg_times, bg_ypos)
        df_bg_IC['timeIC'] = bg_times
        df_bg_IC['locIC'] = bg_locs
        df_bg_C['timeC'] = stim_times
        df_bg_C['locC'] = stim_locs
    print("Flags: ", bsflag1)
    print("Flags: ", blflag1)
    #runtime = bsflag1[-1][-1] + 60
    for i in range(int(len(glut_stim_times))):
        stim_type.append({
            'stim': args.typeStim })
    bg_end = bg_times[-1]
else:
    runtime = 200
trial = args.trial
if chemTurnOff:
    tag_string = 'CytConc_' + str(args.cyt_conc) + '_freq_i_' + str(args.freq_interleaved) + '_spine_create_delay_' + str(args.spine_create_delay) + '_spine_spacing_' + str(args.dummy_spine_spacing)
else:
    tag_string = 'NstimL_' + str(args.NstimL) + '_kcat_enz_tiam_' + str(kcat_enz_tiam) + '_spine_spacing_' + str(args.dummy_spine_spacing) + '_spacing_' + str(args.spacing) + '_bfreq_' + str(args.bfreq) + '_coactive_' + str(args.coactive) + '_mixed_' + str(args.mixed) + '_trial_' + str(args.trial)
df_bg_C.to_csv('./' + tag_string + 'stim.csv')
df_bg_IC.to_csv('./' + tag_string + 'bg.csv')

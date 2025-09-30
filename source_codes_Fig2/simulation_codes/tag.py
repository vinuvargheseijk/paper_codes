from functools import reduce
import sys
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('--kagg', type = float)
parser.add_argument('--ns', type = float)
parser.add_argument('--mu0', type = float)
parser.add_argument('--dt', type = float)
parser.add_argument('--spacing', type = float)
parser.add_argument('--phi_entire', type = float)
parser.add_argument('--cytD', type = float)
args = parser.parse_args()
Length = 20
Kb = 40.0
runtime = 300
bigDia = 1.0
smallDia = 1.0
comptLen = 0.1
#peclet = args.pe


def calcTaper():
  taper = (bigDia - smallDia) / (Length / comptLen)
  return taper


mem_diffConst = 0.0
cyt_diffConst =  args.cytD * 1e-12

#cyt_motor_rate = cyt_diffConst * peclet / (bigDia * 1e-6)
#print("motor rate: ", cyt_motor_rate)
num_spines = int(args.ns)   #num_spines are num_clusters here
cluster_size = 1
cluster_spacing = 1.0 * 1e-6
scalF = 1.0
initial_phitot = 0.5 * args.phi_entire #runs with spine initiation
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
phi_entire = args.phi_entire
dt = args.dt
for nl in range(0, num_spines):
    for cl in range(cluster_size):
      stim_locations.append( start_spine + nl * spacing + cl * cluster_spacing )
      stim_times.append( dt )    
print("STIM LOCATIONS INITIAL: ", stim_locations)      
direction = "forward"    
if direction == "reverse":    
   stim_locations.reverse()
perturb_num = 200
listToStr = reduce(lambda a, b : str(a)+str(b),stim_locations)
print(listToStr)
if len(stim_locations) == 1:
    listToStr = 'One'
tag_string = "phi_entire_" + str(phi_entire) + "_k_agg_" + str(k_agg) + "_ns_" + str(num_spines) + "_mu0_" + str(mu0_exp) + "_spacing_" + str(spacing) + "_dt_" + str(dt) + "_cytD_" + str(cyt_diffConst)
print("Stim locations: ", stim_locations)

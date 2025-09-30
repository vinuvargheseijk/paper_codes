import pandas as pd

#mem_diffConst = 0.01e-12
cyt_diffConst = 1e-12
mu0_exp = -3.5
phi_entire = 9000
num_spines = 2
scalF = 0.65
initial_phitot = 1000
#spacing = 1.0e-6
Length = 10
k_agg = -22

spacings = [0.5e-6, 1.0e-6]
memDs = [0.001e-12, 0.01e-12, 0.1e-12]


for s in range(len(spacings)):
    for d in range(len(memDs)):
        tag_string = "phiE_" + str(phi_entire) + "_ns_" + str(num_spines) + "_mD_" + str(memDs[d]) + "_mu0_" + str(mu0_exp) + "_spacing_" + str(spacings[s]) + "_cD_" + str(cyt_diffConst) + "_k_agg_" + str(k_agg)
        phi = pd.read_csv(tag_string + "phi.csv")
        theta = pd.read_csv(tag_string + "theta.csv")
        print("Final phi: ", phi.iloc[-1])
        print("Num spines: ", len(phi.iloc[-1]))


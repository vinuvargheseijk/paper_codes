import numpy as np
import pandas as pd
import kf_calc



def extract_minEm(phi_entire, file):
    Length = 20
    dendDia = 1e-6
    thickness = 5e-9
    V_cyt = np.pi * ( 0.5 * dendDia )**2 * Length * 1e-6
    Kb = 40.0
    Na = 6.022140e23
    read = file
    df = pd.read_csv(read)
    he = list(df["highEnergy"])
    pe = list(df["phitot"])
    min_index = he.index(np.nanmin(he))
    return df.iloc[he.index(np.nanmin(he))]["highEnergy"], pe[min_index]

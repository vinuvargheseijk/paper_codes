import pandas as pd
import numpy as np



def extract_data(mu0, k_agg1, N, conc):
    df = pd.read_csv("./analytical_N1_N2/N" + str(N) + "_" + "probability_" + str(k_agg1) + "_" + str(conc) + ".csv")
    HRM = float(df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["HRM"])
    HTHETA = float(df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["HTHETA"])
    HRP = float(df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["HRP"])
    LRM = float(df[df["lowEnergy"] == np.nanmin(df["lowEnergy"])]["LRM"])
    LTHETA = float(df[df["lowEnergy"] == np.nanmin(df["lowEnergy"])]["LTHETA"])
    LRP = float(df[df["lowEnergy"] == np.nanmin(df["lowEnergy"])]["LRP"])
    return N * float(df[df["highEnergy"] == np.nanmin(df["highEnergy"])]["highEnergy"]), HRM, HTHETA, HRP, LRM, LTHETA, LRP
 

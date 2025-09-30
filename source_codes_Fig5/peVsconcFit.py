import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("./K-45.0_dE.csv")

phi_entire_list = df["phi_entire"]
dE_list = df["dE"]

Length = 20e-6
diffL = 20e-9

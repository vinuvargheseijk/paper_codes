import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set_context("poster")


filename = "NstimL_10.0_kcat_enz_tiam_50.0_spine_spacing_30.0_spacing_1.0_bfreq_10.0_coactive_1.0_mixed_True_trial_1saved_conc.csv"

df = pd.read_csv(filename)

for i in range(1, 11):
    plt.plot(df["conc" + str(i)])

plt.show()

import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("poster")

kcat = 50.0
bfreq = 10.0
trial = 0

def make_tag(kcat, bfreq, trial):
   f = "NstimL_10.0_kcat_enz_tiam_" + str(kcat) + "_spine_spacing_30.0_spacing_1.0_bfreq_" + str(bfreq) + "_coactive_1.0_mixed_True_trial_" + str(trial)
   return f

tag_string = make_tag(kcat, bfreq, trial)
df = pd.read_csv(tag_string + "saved_conc.csv")
df_marker = pd.read_csv(tag_string + "saved_marker.csv")


fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)
spines = list(df.columns[1:])
colors = ["r", "b", "g", "y", "c", "m"]
axes = [ax1, ax2, ax3, ax4, ax5, ax6]
count = 0
for sp in spines:
    if len(axes) > 0:
      time = np.linspace(0, len(df[sp]) * 0.2, len(df[sp]))
      axes[0].plot(time, df[sp], color = colors[0])
      axes[0].plot(time, df_marker[sp], color = colors[0], linestyle = "-.")
      axes[0].plot(time, [0.4] * len(df[sp]), 'k')
      axes[0].set_title("Spine : " + str(count))
      colors.pop(0)
      axes.pop(0)
      count = count + 1

ax1.set_ylabel("Conc $\mu M$")
ax5.set_xlabel("time (s)")

plt.show()



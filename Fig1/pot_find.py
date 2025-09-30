import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

KBT = 1.38e-23 * 300
file = input("filename")
df = pd.read_csv(file)

CP = list(np.asarray(df["minCP"]) / KBT)
MP = list(np.asarray(df["minMP"]) / KBT)

total = np.asarray(CP) + np.asarray(MP)

plt.plot(df["phitot"], total)
plt.plot(df["phitot"], CP, label = "CP")
plt.plot(df["phitot"], MP, label = "MP")
plt.legend()
plt.figure(101)

plt.plot(df["phitot"], np.asarray(total) / np.asarray(df["phitot"]))

plt.show()


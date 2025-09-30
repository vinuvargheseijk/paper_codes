import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt_times = [1, 1.5, 2, 2.5, 3.5]

def voxelwise_plot(sn):
        temp = []
        for i in plt_times:
            try:
                df = pd.read_csv("./Rd" + str(sn) + str(i) + ".csv")
                plt.plot(df["Rd"])
            except:
                print("File not found")
        plt.ylabel("$R^{d}$")
        plt.xlabel("Voxels")
        plt.show()

import pandas as pd
import matplotlib.pyplot as plt



def get_data(axis, kcat, bfreq, stim_time):
    dfS = pd.read_csv("./time_density_" + str(kcat) +"_S" + str(bfreq) +".csv")
    dfAS = pd.read_csv("./time_density_" + str(kcat) +"_AS" + str(2 * bfreq) +".csv")
    axis.plot(dfS["time"], dfS["density"])
    axis.plot(dfAS["time"], dfAS["density"])
    axis.set_xlim(0, stim_time)
    return dfS, dfAS

fig, ax = plt.subplots()
dfS, dfAS = get_data(ax, 200.0, 10.0, 80.0)
plt.show()

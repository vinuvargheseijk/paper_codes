import numpy as np
import random
import pandas as pd

def poisson_bg(bg_sources, num_sources, freq, coherent):
    initT = 0
    original_sources = bg_sources.copy()
    total_time = initT + num_sources * (1/freq) + 1
    current_total = total_time
    poisson_dt = 1 / freq
    num_samples = num_sources
    df = pd.DataFrame()
    if coherent == False:
      try:
        ncoherent_df = pd.read_csv("./" + str(freq) + "_ncoherent.csv")
        columns = ncoherent_df.columns
        last_column = columns[-1]
        print("dc existing: ", last_column)
        exist = True
      except:
        dcount = 0
        exist = False
      while current_total <= total_time:
         exp_times = np.random.exponential(poisson_dt, 1000)
         print("Total time: ", total_time)
         bg_times = []
         rand_locs = []
         bg_sources = original_sources.copy()
         for et in exp_times[0 : num_samples]:
            temp_locs_index = random.randint(0, len(bg_sources) - 1)
            bg_times.append(initT + et)
            rand_locs.append(bg_sources[temp_locs_index])
            bg_sources.remove(bg_sources[temp_locs_index])
            initT = bg_times[-1]
         current_total = bg_times[-1]
         print("Current duration: ", current_total)
      if exist:
         df[str(int(last_column) + 1)] = bg_times
         df[str( int(last_column) + 1 )] = df[str(int(last_column) + 1)].reset_index(drop=True)
         new_data = pd.concat([ncoherent_df, df], axis = 1)
         new_data.reset_index(drop = True)
         #ncoherent_df[str(int(last_column) + 1)] = bg_times
         new_data.to_csv("./" + str(freq) + "_ncoherent.csv")
      else:
         df["0"] = bg_times
         df.to_csv("./" + str(freq) + "_ncoherent.csv")
    if coherent == True:
      try:
        coherent_df = pd.read_csv("./" + str(freq) + "_coherent.csv")
        columns = coherent_df.columns
        last_column = columns[-1]
        print("dc existing: ", last_column)
        exist = True
      except:
        dcount = 0
        exist = False
      while current_total >= total_time:
         exp_times = np.random.exponential(poisson_dt, 1000)
         bg_times = []
         rand_locs = []
         bg_sources = original_sources.copy()
         for et in exp_times[0 : num_samples]:
            temp_locs_index = random.randint(0, len(bg_sources) - 1)
            bg_times.append(initT + et)
            rand_locs.append(bg_sources[temp_locs_index])
            bg_sources.remove(bg_sources[temp_locs_index])
            initT = bg_times[-1]
         current_total = bg_times[-1]
      if exist:
         df[str(int(last_column) + 1)] = bg_times
         df[str( int(last_column) + 1 )] = df[str(int(last_column) + 1)].reset_index(drop=True)
         new_data = pd.concat([coherent_df, df], axis = 1)
         new_data.reset_index(drop = True)
         #ncoherent_df[str(int(last_column) + 1)] = bg_times
         new_data.to_csv("./" + str(freq) + "_coherent.csv")
      else:
         df["0"] = bg_times
         df.to_csv("./" + str(freq) + "_coherent.csv")
    return (bg_times, rand_locs)

num_stim = 5
locations = list(np.linspace(0,10, num_stim))

print(poisson_bg(locations, num_stim, 20.0, False))
print(poisson_bg(locations, num_stim, 20.0, True))

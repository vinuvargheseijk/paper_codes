import numpy as np
import matplotlib.pyplot as plt


def vol_cap(rd, rs, theta):
    a = (rd + rs) * np.sin(theta)
    h = rd * (1 - np.cos(theta))
    vol = (1/6.0) * np.pi * h * (3 * a**2 + h**2)
    return vol


def shank_reac(phitot, shank_conc, rd, rs, theta, kf, kb):
   dt = 0.01
   a_list = []
   b_list = []
   ab_list = []
   volume = vol_cap(rd, rs, theta)
   a0 = 1e3 * phitot / (6.0223e23 * volume)
   print("Conc membrane (uM): ", a0)
   print("Volume: ", volume)
   b0 = shank_conc
   ab0 = 0
   for i in range(10000):
       ab_update = ab0 + dt * (kf * a0 * b0 - kb * ab0)
       b_update = b0 + dt * (-kf * a0 * b0 + kb * ab0)
       a_update = a0 + dt * (-kf * a0 * b0 + kb * ab0)
       a_list.append(a_update)
       b_list.append(b_update)
       ab_list.append(ab_update)
       ab0 = ab_update
       a0 = a_update
       b0 = b_update
   return a_list[-1], b_list[-1], ab_list[-1]  


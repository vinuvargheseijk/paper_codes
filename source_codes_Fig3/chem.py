import numpy as np
import matplotlib.pyplot as plt


a0 = 0.1
b0 = 0.5
ab0 = 0
kf = 0.5
kb = 0.1
dt = 0.01

a_list = []
b_list = []
ab_list = []

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
    
plt.plot(ab_list, label = "ab")    
plt.plot(b_list, label = "b")    
plt.plot(a_list, label = "a")    
plt.plot( (kb/kf) * np.asarray(ab_list) / np.asarray(a_list) - np.asarray(b_list), linestyle = None, marker = "*" )

plt.legend()

plt.show()

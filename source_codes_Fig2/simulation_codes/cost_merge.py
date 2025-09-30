import numpy as np
import matplotlib.pyplot as plt
import energy_process as EP
from scipy import optimize
from scipy.optimize import LinearConstraint, NonlinearConstraint



thickness = EP.thickness
fs = EP.fs
phisat =  EP.phisat
KBT = EP.KBT
Na = EP.Na
RT = EP.RT
RT_kcal_mol = EP.RT_kcal_mol #kcal/mol
k = EP.k
khat = EP.khat
kentropy = EP.KBT
Cp = EP.Cp
dendDia = EP.dendDia
Rc = EP.Rc
k_agg = -50.0
mu0_exp = -4.0
phitot = 1500.0
phi_entire = 4000.0
conc = phi_entire-phitot

def calc_saddle_area(alpha, rs, rm, theta, frac):
    r0 = (rm + rs) * np.sin(alpha)
    area_saddle = 2 * np.pi * rs * (r0 * theta + rs * (np.cos(alpha) - 1))
    return frac * area_saddle * fs

def pull_cost(alpha, rs, rm, theta, frac):
    yr = rs * (1 - np.cos(alpha))
    area_saddle = calc_saddle_area(alpha, rs, rm, theta, frac)
    return area_saddle * ( (yr - Rc)**2 / Rc**2 ) * fs


def tension_relax(rs, rm, theta, alpha, frac):
    frac_area = calc_saddle_area(alpha, rs, rm, theta, frac)
    return frac_area * fs / 2



def calc_phi(xp): 
    area = EP.calc_dome_area( (xp[0]/xp[1]) * 1e-6, xp[1])
    phi = phitot / (area * EP.phisat)
    return phi

linear_constraint1 = LinearConstraint( [[1, 0, 0]], [0.001], [5.0] )
linear_constraint2 = LinearConstraint( [[0, 1, 0]], [0.01], [1.5] )
linear_constraint3 = LinearConstraint( [[0, 0, 1]], [0.001], [5.0] )
nonlinear_constraint = NonlinearConstraint(calc_phi, 0.2, 1.0)

sample = [0.5, 0.5, 0.5]
ret = optimize.minimize( EP.total_energy, sample, method = 'trust-constr', args = (phitot, conc, k_agg, mu0_exp), constraints = [linear_constraint1, linear_constraint2, linear_constraint3, nonlinear_constraint])



print(ret)
rm = ret.x[0] / ret.x[1]
rm = rm * 1e-6
theta = ret.x[1]
rs = ret.x[2] * 1e-6

alpha = np.linspace(0, theta, 20)

trelax = []
pcost = []
frac = 0.5

for a in alpha:
    trelax.append(tension_relax(rs, rm, theta, a, frac) / KBT)
    pcost.append(pull_cost(a, rs, rm, theta, frac)/ KBT)


plt.plot(alpha, trelax, label = "Relax")
plt.plot(alpha, pcost, label = "cost")
plt.legend()
plt.show()

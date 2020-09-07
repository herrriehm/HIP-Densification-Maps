"""
    HIP 20.0 - calculate densification maps for hot isostatic pressing.
    Copyright (C) 2020 Sebastian Riehm

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import matplotlib.pyplot as plt
# import numpy as np
from scipy.integrate import solve_ivp

from fun_aux import *
from fun_dens import *

STEPS = 50

TMIN = 1200 + 273.0
TMAX = 1200 + 273.0

Pmax = 100.0E6  # Pa

TIMES = np.array([15, 30, 60, 120, 240])  # min
TIME = np.max(TIMES)

TArray = np.linspace(TMIN, TMAX, num=STEPS)
PArray = np.logspace(-2, 1, num=STEPS) * material.SIGMAY
# PArray = np.linspace(150.0E6, 150.0E6, num=STEPS)

# calculate initial density Dy from instant yielding
DYield = np.zeros(STEPS)

for i in np.arange(STEPS):
    DYield[i] = d_yield1(TArray[i], PArray[i])
    if DYield[i] >= constants.DBREAK:
        DYield[i] = d_yield2(TArray[i], PArray[i])

t_eval = np.arange(0, TIME * 60 + 1, 300)

DTotal = np.zeros([STEPS, 5])

for i in np.arange(STEPS):
    DStart = DYield[i]
    P = PArray[i]
    T = TArray[i]
    if material.G < 2 * x(DStart):
        sol = solve_ivp(d_total_with_nhc, [0, TIME * 60], [DStart],
                        args=[T, P], t_eval=t_eval)
        DStep = sol.y[0][[3, 6, 12, 24, 48]]
    else:
        sol = solve_ivp(d_total_without_nhc, [0, TIME * 60], [DStart],
                        args=[T, P], t_eval=t_eval)
        DStep = sol.y[0][[3, 6, 12, 24, 48]]

    DTotal[i] = DStep

fig, ax = plt.subplots()
ax.semilogx(PArray, DYield)
for i in np.arange(np.size(TIMES)):
    ax.semilogx(PArray, DTotal[:, i])
# ax.semilogx([np.amin(PArray), np.amax(PArray)], [constants.D0, constants.D0])
plt.xlim(np.amin(PArray), np.amax(PArray))
plt.ylim(0.6, 1)
plt.xlabel('Druck / Pa')
plt.ylabel('relative Dichte / –')
plt.show()

# plt.plot(TArray, DYield)
# plt.plot(TArray, DTotal)
# plt.xlim(TArray[0], TArray[-1])
# plt.ylim(0.985, 1)
# plt.xlabel('Temperatur / K')
# plt.ylabel('relative Dichte / –')
# plt.show()

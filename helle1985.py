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
from scipy.integrate import solve_ivp

from fun_aux import *
from fun_dens import *

TMIN = 1673
TMAX = 1723

Pmax = 100.0E6  # Pa

# times for contour lines
TIMES = np.array([15, 30, 60, 120, 300])  # min

# total time of HIP cycle to be simulated
TIME = np.max(TIMES)

# times at which to store the computed solution of rate equations
# every minute seems to be okay
t_eval = np.arange(0, TIME * 60 + 1, 60)

# number of steps 'on the x axis' to calculate the densification
STEPS = np.size(t_eval)

TArray = np.linspace(TMIN, TMAX, num=STEPS)

def TArrayRedouani2019():
    if t<60:
        pass
    else:
        pass
    return


# Data from Arzt 1983 Figure 3
# TArray = np.array(
#     [21.20777012, 61.84604206, 102.4836087, 143.1264648, 183.7636789,
#      224.3986009, 265.0363439, 305.6707369, 346.3046011, 386.9358205,
#      427.5702136, 468.195791, 508.8062055, 549.435726, 590.0264175,
#      630.6455596, 671.2573846, 711.8611634, 752.4650383, 793.0685286,
#      834.6345424, 942.9579064, 996.066297, 1036.65573, 1077.247456,
#      1117.832707, 1158.409749, 1198.989661, 1239.56305, 1280.136439,
#      1320.69872, 1361.254478, 1401.807767, 1442.356649, 1482.885607,
#      1523.403457, 1563.915842, 1604.412359, 1644.897415, 1678.21225])

# PArray = np.logspace(-3, 1, num=STEPS) * material.SIGMAY
PArray = np.linspace(Pmax, Pmax, num=STEPS)

# calculate initial density Dy from instant yielding
DYield = np.zeros(STEPS)

# Data from Arzt 1983 Figure 3
# DYield = np.array(
#     [0.667967908, 0.667949179, 0.667949238, 0.667808394, 0.667817846,
#      0.667888356, 0.667883718, 0.667968319, 0.66806701, 0.668236153,
#      0.668320753, 0.668640193, 0.669363554, 0.669577957, 0.670826717,
#      0.671317588, 0.672003376, 0.672903507, 0.673801077, 0.674708893,
#      0.675372978, 0.67864781, 0.679772309, 0.681054586, 0.682275806,
#      0.683669478, 0.685281858, 0.686817762, 0.688527446, 0.69023713,
#      0.692242711, 0.694422072, 0.696667188, 0.699029724, 0.701922995,
#      0.705112163, 0.70844693, 0.712204407, 0.716267175, 0.719273171])

for i in np.arange(STEPS):
    DYield[i] = d_yield1(TArray[i], PArray[i])
    if DYield[i] >= constants.DBREAK:
        DYield[i] = d_yield2(TArray[i], PArray[i])

t_eval = np.arange(0, TIME * 60 + 1, 60)

DTotal = np.zeros([STEPS, 5])
DTime = np.zeros([STEPS])

for i in np.arange(STEPS):
    DStart = DYield[i]
    P = PArray[i]
    T = TArray[i]
    if material.G < 2 * x(DStart):
        sol = solve_ivp(d_total_with_nhc, [0, TIME * 60], [DStart],
                        args=[T, P], t_eval=t_eval)
    else:
        sol = solve_ivp(d_total_without_nhc, [0, TIME * 60], [DStart],
                        args=[T, P], t_eval=t_eval)

    DStep = sol.y[0][[15, 30, 60, 120, 240]]
    DTotal[i] = DStep
    DTime[i] = sol.y[0][i]


# fig, ax = plt.subplots()
# ax.semilogx(PArray, DYield)
# plt.plot(PArray, DTotal)
# # ax.semilogx([np.amin(PArray), np.amax(PArray)], [constants.D0, constants.D0])
# plt.xlim(np.amin(PArray), np.amax(PArray))
# plt.ylim(0.6, 1)
# plt.xlabel('Druck / Pa')
# plt.ylabel('relative Dichte / –')
# plt.show()

# def fitjasper(temperature):
#     return 0.9973 * np.exp(-5.552e-08 * temperature) - 2.106e+04 * np.exp(
#         -0.008148 * temperature)
#
#
# DJasper = fitjasper(TArray)

plt.plot(TArray, DYield)
plt.plot(TArray, DTotal)
plt.xlim(TArray[0], TArray[-1])
plt.ylim(0.7, 1)
plt.xlabel('Temperatur / K')
plt.ylabel('relative Dichte / –')
plt.show()

csvArray = np.concatenate(
    (np.transpose([TArray]), np.transpose([DYield]), DTotal), axis=1)

exportfilename = 'intermediateForDoverT.csv'

np.savetxt(exportfilename, csvArray,
           header='T, Yield, 15 min, 30 min, 60 min, 120 min, 240 min',
           delimiter=',')

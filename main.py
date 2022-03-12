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

from fun_dens import *

# from ashby1987 import *

# import reardon1998
# import svoboda1996
# import schnecke

import cvntable58

# from redouani2012 import *

# TMIN = 25 + 273
# TMAX = 25 + 273
# Pmin = 250E6
# Pmax = 250E6

TArray = cvntable58.TArray
PArray = cvntable58.PArray

# times for contour lines
TIMES = TArray.size-1 #np.array([60])  # min

# total minutes of HIP cycle to be simulated
TIME = np.max(TIMES)

# times at which to store the computed solution of rate equations
# every minute seems to be okay
t_eval = np.arange(0, TIME * 60 + 1, 60)

# number of steps 'on the x axis' to calculate the densification
STEPS = np.size(t_eval)

# TArray = np.linspace(TMIN, TMAX, num=STEPS)
# PArray = np.logspace(-2, 1,
#                      num=STEPS) * fun_aux.sigmay(TMAX)
# PArray = np.linspace(Pmin, Pmax, num=STEPS)

# calculate initial density Dy from instant yielding
DYield = np.zeros(STEPS)

for i in np.arange(STEPS):
    DYield[i] = d_yield1(TArray[i], PArray[i])
    if DYield[i] >= constants.DBREAK:
        DYield[i] = d_yield2(TArray[i], PArray[i])

DTotal = np.ones([STEPS, np.size(TIMES)])
DTime = np.zeros([STEPS])


def fulldense(t, y, temperature, pressure):
    """event function for solve_ivp to check if D is already 1"""
    return y[0] - 1


fulldense.terminal = True

for i in np.arange(STEPS):
    DStart = DYield[i]
    P = PArray[i]
    T = TArray[i]
    sol = solve_ivp(d_total, [0, TIME * 60], [DStart],
                    args=[T, P], t_eval=t_eval, rtol=1e-6, atol=1e-6)

    DStep = sol.y[0][TIMES]
    DTotal[i] = DStep
    DTime[i] = sol.y[0][i]
    # print(DTime[i])

DTotal = np.where(DTotal <= 1, DTotal, 1)
DTime = np.where(DTime <= 1, DTime, 1)


# fig, ax = plt.subplots(figsize=(8, 8))
# ax.semilogx(PArray, DYield)
# # # # plt.plot(TArray, DYield)
# plt.plot(PArray, DTotal)
# # # # # ax.semilogx([np.amin(PArray), np.amax(PArray)], [constants.D0, constants.D0])
# plt.xlim(np.amin(PArray), np.amax(PArray))
# plt.ylim(0.6, 1)
# plt.xlabel('Pressure / Pa')
# plt.ylabel('Relative Density / –')
# plt.show()

# print(material.R)
# def fitjasper(temperature):
#     return 0.9973 * np.exp(-5.552e-08 * temperature) - 2.106e+04 * np.exp(
#         -0.008148 * temperature)
#
#
# DJasper = fitjasper(TArray)

print(DTime[-1])

# fig, ax = plt.subplots()
# # plt.plot(TArray, DYield)
# # plt.plot(TArray, DTotal)
# plt.plot(DTime, label='D')
# # ax.plot(TArray/1000, label='T')
# # ax.plot(PArray/100000000, label='P')
# #
# # legend = ax.legend(loc='upper center', fontsize='x-large')
#
# # plt.xlim(TMIN, TMAX)
# plt.ylim(0.64, 1)
# plt.xlabel('Zeit / min')
# plt.ylabel('relative Dichte / –')
# plt.show()
# print(f'T: {TArray[-1]-273}, P: {PArray[-1]/1000000}, D: {DTime[-1]}')
#
# csvArray = np.concatenate(
#     (np.transpose([PArray]), np.transpose([DYield]), DTotal), axis=1)

# exportfilename = "100.csv"
#
# np.savetxt(material.exportfilename, DTotal*100,
# #            # header='P, Yield, 15 min, 30 min, 60 min, 120 min, 240 min',
#            delimiter=',')

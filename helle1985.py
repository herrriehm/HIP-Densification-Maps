import matplotlib.pyplot as plt
import numpy as np
import scipy.constants
from scipy.integrate import solve_ivp

# Tool steel
OMEGA = 1.21E-29  # m^3
B = 2.58E-10  # m
TM = 1680  # K
SIGMAY = 200.0E6  # Pa
SIGMAY0 = 9.1E8  # Pa
ALPHA = -4.5
MU0 = 8.1E10  # Pa
BETA = -0.85
D0V = 3.7E-5  # m^2 s^-1
QV = 280.0E3  # J mol^-1
DD0B = 2.0E-13  # m^3 s^-1
QB = 167.0E3  # J mol^-1
n = 7.5
A = 1.5E12

R = 25.0E-6  # m
G = 10.0E-6  # m

Z0 = 7.3
C = 15.5

D0 = 0.64
DBREAK: float = 0.9
STEPS = 15

CONSTR = scipy.constants.R  # J/mol/K
CONSTK = scipy.constants.k  # J/K

TMIN = 1200 + 273.0
TMAX = 1200 + 273.0

Pmax = 100.0E6  # Pa

TIME = 60  # min


def mu(t):
    return MU0 * (1 + (t - 300) / TM * BETA)


def sigmay(t):
    return SIGMAY0 * (1 + (t - 300) / TM * ALPHA)


def dv(t):
    return D0V * np.exp(-QV / (CONSTR * t))


def db(t):
    return DD0B * np.exp(-QB / (CONSTR * t))


def epssig(t):
    return A * B * dv(t) / (CONSTK * t * np.power(mu(t), n - 1))


def x(d):
    return np.power(3, -1 / 2) * np.power((d - D0) / (1 - D0), 1 / 2) * R


def p_eff(d):
    return P * (1 - D0) / (np.power(d, 2) * (d - D0))


def rho(d):
    return R * (d - D0)


def r(d):
    if d > 1:
        d = 1
    return R * np.power((1 - d) / 6, 1 / 3)


def d_yield1(p, t):
    return np.power((1 - D0) * p / (1.3 * sigmay(t)) + np.power(D0, 3), 1 / 3)


def d_yield2(p, t):
    return 1 - np.exp(-3 / 2 * p / sigmay(t))


def d_dot_plc1(d):
    if d > 1:
        d = 1
    return 5.3 * np.power(np.power(d, 2) * D0, 1 / 3) * x(d) / R * epssig(T) * np.power(p_eff(d) / 3, n)


def d_dot_plc2(d):
    if d > 1:
        d = 1
    return 3 / 2 * epssig(T) * d * (1 - d) / np.power(1 - np.power(1 - d, 1 / n), n) * np.power(3 / 2 * P / n, n)


def d_dot_ipb1(d):
    if d > 1:
        d = 1
    return 43 * np.power(1 - D0, 2) / np.power(d - D0, 2) * (db(T) + rho(d) * dv(T)) / (
            CONSTK * T * np.power(R, 3)) * OMEGA * P


def d_dot_ipb2(d):
    if d > 1:
        d = 1
    return 54 * OMEGA * (db(T) + r(d) * dv(T)) / (CONSTK * T * np.power(R, 3)) * 5 * np.power(1 - d, 1 / 2) * P


def d_dot_nhc1(d):
    if d > 1:
        d = 1
    return 24.9 * OMEGA / (CONSTK * T * np.power(G, 2)) * np.power(np.power(d, 2) * D0, 1 / 3) * x(d) / R * (
            dv(T) + np.pi * db(T) / G) * p_eff(d)


def d_dot_nhc2(d):
    if d > 1:
        d = 1
    return 31.5 * OMEGA / (CONSTK * T * np.power(G, 2)) * (1 - d) * (dv(T) + np.pi * db(T) / G) * P


def d_total_without_nhc1(t, d):
    return d_dot_plc1(d) + d_dot_ipb1(d)


def d_total_without_nhc2(t, d):
    return d_dot_plc2(d) + d_dot_ipb2(d)


def d_total_with_nhc1(t, d):
    return d_dot_plc1(d) + d_dot_ipb1(d)  # + d_dot_nhc1(d)


def d_total_with_nhc2(t, d):
    return d_dot_plc2(d) + d_dot_ipb2(d)  # + d_dot_nhc2(d)


TArray = np.linspace(1, TM, num=STEPS)
# PArray = np.logspace(-2, 1, num=STEPS) * SIGMAY
PArray = np.linspace(100.0E6, 100.0E6, num=STEPS)

# calculate initial density Dy from instant yielding
DYield = np.zeros(STEPS)

for i in np.arange(STEPS):
    DYield[i] = d_yield1(PArray[i], TArray[i])
    if DYield[i] >= DBREAK:
        DYield[i] = d_yield2(PArray[i], TArray[i])

t_eval = np.arange(0, TIME * 60 + 1, 300)

print(DYield)

# DTotal = np.zeros(STEPS)
# # Dipb = np.zeros(STEPS)
# # # Dnhc=np.zeros(STEPS)
# # Dgesamt = np.zeros(STEPS)
# #
# # # zx = np.zeros(STEPS)
# #
# # T = 1200 + 273.0
# #
# for i in np.arange(STEPS):
#     # print(i)
#     DStart = DYield[i]
#     P = PArray[i]
#     T = TArray[i]
#     if DStart < DBREAK:
#         if G < 2 * x(DStart):
#             DStep = solve_ivp(d_total_with_nhc1, [0, TIME * 60], [DStart], t_eval=t_eval).y[0][-1]
#         else:
#             DStep = solve_ivp(d_total_without_nhc1, [0, TIME * 60], [DStart], t_eval=t_eval).y[0][-1]
# #             # Dipbstep = solve_ivp(Dipb1, [0, TIME*60], [Dstart]).y[0][-1]
# #             # if G > 2*x(Dstart):
# #             #    Dnhcstep = solve_ivp(Dnhc1, [0, TIME*60], [Dstart]).y[0][-1]
# #             #    print('G>2x in 1!')
# #             # else:
# #             #    Dnhcstep = Dy[i]
# #             #    print('G<2x in 1!')
#     else:
#         if G < 2 * x(DStart):
#             DStep = solve_ivp(d_total_with_nhc2, [0, TIME * 60], [DStart], t_eval=t_eval).y[0][-1]
#         else:
#             DStep = solve_ivp(d_total_without_nhc2, [0, TIME * 60], [DStart], t_eval=t_eval).y[0][-1]
# #
# #             # Dplcstep = solve_ivp(Dplc2, [0, TIME*60], [Dstart]).y[0][-1]
# #             # Dipbstep = solve_ivp(Dipb2, [0, TIME*60], [Dstart]).y[0][-1]
# #             # if G > 2*x(Dstart):
# #             #    Dnhcstep = solve_ivp(Dnhc2, [0, TIME*60], [Dstart]).y[0][-1]
# #             #    print('G>2x in 2!')
# #             # else:
# #             #    Dnhcstep = Dy[i]
# #             #    print('G<2x in 2!')
#     if DStep > 1:
#         DStep = 1.0
# #     # print(Dplcstep)
#     DTotal[i] = DStep
# #
# # #    zx[i] = 2*x(Dstart)
# #
# # # if Dplcstep > 1: Dplcstep = 1
# # # if Dipbstep > 1: Dipbstep = 1
# # #    if Dnhcstep > 1: Dnhcstep = 1
# # # Dplc[i] = Dplcstep - Dy[i]
# # # Dipb[i] = Dipbstep - Dy[i]
# # #    Dnhc[i] = Dnhcstep - Dy[i]
#
# # fig, ax = plt.subplots(constrained_layout=True)
# # ax.semilogx(PArray, DYield)
# # ax.semilogx(PArray, DTotal)
# # plt.xlim(np.amin(PArray), np.amax(PArray))
# # plt.ylim(0.6, 1)
# # plt.xlabel('Druck / Pa')
# # plt.ylabel('relative Dichte / –')
# # plt.show()
#

x = np.arange(1, 1680)
y = mu(x)

plt.plot(x, y)
# plt.xlim(1, TM)
# plt.ylim(0.6, 1)
plt.xlabel('Temperatur / K')
# plt.ylabel('relative Dichte / –')
plt.show()

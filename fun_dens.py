import numpy as np

import constants
import fun_aux
import material


# def d_yield1(temperature, pressure):
#     return np.power((1 - constants.D0) * pressure / (1.3 * fun_aux.sigmay(temperature)) +
#                     np.power(constants.D0, 3), 1 / 3)
#
#
# def d_yield2(temperature, pressure):
#     return 1 - np.exp(-3 / 2 * pressure / fun_aux.sigmay(temperature))


def d_yield1(temperature, pressure):
    return np.power((1 - constants.D0) * pressure / (1.3 * material.SIGMAY) +
                    np.power(constants.D0, 3), 1 / 3)


def d_yield2(temperature, pressure):
    return 1 - np.exp(-3 / 2 * pressure / material.SIGMAY)


def d_dot_plc1(d, temperature, pressure):
    return 5.3 * np.power(np.power(d, 2) * constants.D0, 1 / 3) * \
           fun_aux.x(d) / material.R * fun_aux.epssig(temperature) * \
           np.power(fun_aux.p_eff(d, pressure) / 3, material.n)


def d_dot_plc2(d, temperature, pressure):
    return 3 / 2 * fun_aux.epssig(temperature) * d * (1 - d) / \
           np.power(1 - np.power(1 - d, 1 / material.n), material.n) * \
           np.power(3 / 2 * pressure / material.n, material.n)


def d_dot_ipb1(d, temperature, pressure):
    return 43 * np.power(1 - constants.D0, 2) / np.power(d - constants.D0, 2) \
           * (fun_aux.db(temperature) + fun_aux.rho(d) *
              fun_aux.dv(temperature)) / (constants.BOLTZMANN * temperature *
                                          np.power(material.R, 3)) * \
           material.OMEGA * pressure


def d_dot_ipb2(d, temperature, pressure):
    return 54 * material.OMEGA * (fun_aux.db(temperature) + fun_aux.r(d) *
                                  fun_aux.dv(temperature)) / (
                   constants.BOLTZMANN * temperature *
                   np.power(material.R, 3)) * 5 * np.power(1 - d,
                                                           1 / 2) * pressure


def d_dot_nhc1(d, temperature, pressure):
    return 24.9 * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.G, 2)) * \
           np.power(np.power(d, 2) * constants.D0, 1 / 3) * fun_aux.x(
        d) / material.R * \
           (fun_aux.dv(temperature) + np.pi * fun_aux.db(
               temperature) / material.G) * fun_aux.p_eff(d, pressure)


def d_dot_nhc2(d, temperature, pressure):
    return 31.5 * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.G,
                                                         2)) * (1 - d) * \
           (fun_aux.dv(temperature) + np.pi * fun_aux.db(
               temperature) / material.G) * pressure


def d_total_without_nhc1(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # print(f'1  D: {d}, T: {temperature}, P: {pressure}')
    return d_dot_plc1(d, temperature, pressure) + \
        d_dot_ipb1(d, temperature, pressure)


def d_total_without_nhc2(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # print(f'2  D: {d}, T: {temperature}, P: {pressure}')
    return d_dot_plc2(d, temperature, pressure) + \
        d_dot_ipb2(d, temperature, pressure)


def d_total_with_nhc1(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # print(f'3  D: {d}, T: {temperature}, P: {pressure}')
    return d_dot_plc1(d, temperature, pressure) + \
        d_dot_ipb1(d, temperature, pressure) + \
        d_dot_nhc1(d, temperature, pressure)


def d_total_with_nhc2(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # print(f'4  D: {d}, T: {temperature}, P: {pressure}')
    return d_dot_plc2(d, temperature, pressure) + \
        d_dot_ipb2(d, temperature, pressure) + \
        d_dot_nhc2(d, temperature, pressure)


def smoothing(d):
    """
    smoothing function between stages 1 and 2 according to
    Messung und Simulation der Verdichtungskinetik pulvermetallurgischer
    Hochtemperaturwerkstoffe beim heißisostatischen Pressen by Martin
    Dietze (1991)
    """

    return 0.5 * (1 + (np.exp(75 * (d - 0.9)) - 1) / (
            np.exp(75 * (d - 0.9)) + 1))


def d_total_without_nhc(t, d, temperature, pressure):
    return (1 - smoothing(d)) * d_total_without_nhc1(d, temperature,
                                                     pressure) + \
           smoothing(d) * d_total_without_nhc2(d, temperature, pressure)


def d_total_with_nhc(t, d, temperature, pressure):
    return (1 - smoothing(d)) * d_total_with_nhc1(d, temperature, pressure) + \
           smoothing(d) * d_total_with_nhc2(d, temperature, pressure)

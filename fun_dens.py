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

import numpy as np

import constants
import fun_aux
import material


# def d_yield1(temperature, pressure):
#     return np.power(
#         (1 - constants.D0) * pressure / (1.3 * fun_aux.sigmay(temperature)) +
#         np.power(constants.D0, 3), 1 / 3)
#
#
# def d_yield2(temperature, pressure):
#     result = 1 - np.exp(-3 / 2 * pressure / fun_aux.sigmay(temperature))
#     if result < constants.D0:
#         return constants.D0
#     else:
#         return result

# def d_yield1(temperature, pressure):
#     """Dyield1 from Redouani 2019, eq. 8"""
#     return constants.D0 / 2 * (1 + np.power(
#         1 + 2 / (3 * constants.D0) * pressure / fun_aux.sigmay(temperature),
#         1 / 2))
#
#
# def d_yield2(temperature, pressure):
#     """Dyield2 from Redouani 2019, eq. 10"""
#     return 1 - np.exp(-3 / 2 * pressure / fun_aux.sigmay(temperature))


def d_yield1(temperature, pressure):
    return np.power((1 - constants.D0) * pressure / (1.3 * material.SIGMAY0) +
                    np.power(constants.D0, 3), 1 / 3)


def d_yield2(temperature, pressure):
    result = 1 - np.exp(-3 / 2 * pressure / material.SIGMAY0)
    if result < constants.D0:
        return constants.D0
    else:
        return result


def d_dot_plc1(d, temperature, pressure):
    return 5.3 * np.power(np.power(d, 2) * constants.D0, 1 / 3) * \
           fun_aux.x(d) / material.R * fun_aux.epssig(temperature) * \
           np.power(fun_aux.p_eff(d, pressure) / 3, material.n)


# def d_dot_plc1(d, temperature, pressure):
#     """modification from McCoy"""
#     return 5.3 * d * np.power(d / constants.D0, 1 / 3) * \
#            fun_aux.x(d) / material.R * fun_aux.epssig(temperature) * \
#            np.power(fun_aux.p_eff(d, pressure) / 3, material.n)


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

#
def d_dot_nhc1(d, temperature, pressure):
    if material.G < 2 * fun_aux.x(d):
        return 24.9 * material.OMEGA / (
                constants.BOLTZMANN * temperature * np.power(material.G, 2)) * \
               np.power(np.power(d, 2) * constants.D0, 1 / 3) * fun_aux.x(
            d) / material.R * \
               (fun_aux.dv(temperature) + np.pi * fun_aux.db(
                   temperature) / material.G) * fun_aux.p_eff(d, pressure)
    else:
        return 0


# def d_dot_nhc1(d, temperature, pressure):
#     """modification from McCoy"""
#     return 24.9 * material.OMEGA / (constants.BOLTZMANN * temperature
#            * np.power(material.G, 2)) * d * np.power(d / constants.D0, 1 / 3) \
#            * fun_aux.x(d) / material.R * (fun_aux.dv(temperature) + np.pi
#            * fun_aux.db(temperature) / material.G) * fun_aux.p_eff(d, pressure)


def d_dot_nhc2(d, temperature, pressure):
    if material.G < 2 * fun_aux.x(d):
        return 31.5 * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.G,
                                                         2)) * (1 - d) * \
           (fun_aux.dv(temperature) + np.pi * fun_aux.db(
               temperature) / material.G) * pressure
    else:
        return 0


# def d_total_without_nhc1(d, temperature, pressure):
#     if d > 1.0:
#         d = 1.0
#     # if d < constants.D0:
#     #     d = constants.D0+1.0E-6
#     return d_dot_plc1(d, temperature, pressure) + \
#            d_dot_ipb1(d, temperature, pressure)
#
#
# def d_total_without_nhc2(d, temperature, pressure):
#     if d > 1.0:
#         d = 1.0
#     # if d < constants.D0:
#     #     d = constants.D0+1.0E-6
#     return d_dot_plc2(d, temperature, pressure) + \
#            d_dot_ipb2(d, temperature, pressure)


def d_total_1(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
    return d_dot_plc1(d, temperature, pressure) + \
           d_dot_ipb1(d, temperature, pressure) + \
           d_dot_nhc1(d, temperature, pressure)


def d_total_2(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
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
    if d > 1.0:
        d = 1.0

    return 0.5 * (1 + (np.exp(75 * (d - constants.DBREAK)) - 1) / (
            np.exp(75 * (d - constants.DBREAK)) + 1))


# def d_total_without_nhc(t, d, temperature, pressure):
#     return (1 - smoothing(d)) * d_total_without_nhc1(d, temperature,
#                                                      pressure) + \
#            smoothing(d) * d_total_without_nhc2(d, temperature, pressure)


def d_total(t, d, temperature, pressure):
    return (1 - smoothing(d)) * d_total_1(d, temperature, pressure) + \
           smoothing(d) * d_total_2(d, temperature, pressure)

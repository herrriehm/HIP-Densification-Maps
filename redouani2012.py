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


def g(d):
    return (np.power(d / constants.D0, 1 / 3) - 1) * (
            2 * constants.Z0 + constants.C * (
            np.power(d / constants.D0, 1 / 3) - 1))


def d_yield1(temperature, pressure):
    """Redouani 2019"""
    result = constants.D0 / 2 * (1 + np.power(
        1 + 2 / (3 * constants.D0) * pressure / fun_aux.sigmay(temperature),
        1 / 2))
    if result < constants.D0:
        return constants.D0
    else:
        return result


def d_yield2(temperature, pressure):
    """Redouani 2019"""
    result = 1 - np.exp(-3 / 2 * pressure / fun_aux.sigmay(temperature))
    if result < constants.D0:
        return constants.D0
    else:
        return result

# def d_yield1(temperature, pressure):
#     """Redouani 2012"""
#     result = constants.D0 / 2 * (1 + np.power(
#         1 + 2 / (3 * constants.D0) * pressure / fun_aux.sigmay(temperature),
#         1 / 2))
#     if result < constants.D0:
#         return constants.D0
#     else:
#         return result
#
#
# def d_yield2(temperature, pressure):
#     """Redouani 2012"""
#     result = constants.D0 / 2 * (1 + np.power(
#         1 + 2 / (3 * constants.D0) * pressure / fun_aux.sigmay(temperature),
#         1 / 2))
#     if result < constants.D0:
#         return constants.D0
#     else:
#         return result


def d_dot_plc1(d, temperature, pressure):
    result = 2 * fun_aux.A(temperature) * np.power(constants.D0 / d,
                                                   1 / 3) * np.power(
        d * (d - constants.D0) / constants.D0, 1 - material.n) * np.power(
        pressure / 6, material.n)
    if result > 1:
        result = 1
    return result


def d_dot_plc2(d, temperature, pressure):
    return 3 / 2 * fun_aux.A(temperature) * d * (1 - d) / np.power(
        1 - np.power(1 - d, 1 / material.n), material.n) * np.power(
        3 / 2 * pressure / material.n, material.n)


def d_dot_gbd1(d, temperature, pressure):
    return 72 * fun_aux.db(temperature) * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.R,
                                                         3)) * np.power(
        constants.D0, 5 / 3) * np.power(d, 1 / 3) / (
                   (d - constants.D0) * g(d)) * pressure


def d_dot_gbd2(d, temperature, pressure):
    return 27 * fun_aux.db(temperature) * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.R,
                                                         3)) * np.power(d,
                                                                        2) / (
                   1 - 5 * np.power((1 - d) / (5 * d), 2 / 3)) * pressure


def d_dot_ld1(d, temperature, pressure):
    return 59 * fun_aux.dv(temperature) * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.R,
                                                         2)) * np.power(
        constants.D0, 7 / 6) * np.power(d, 1 / 3) / (
                   np.power(d - constants.D0, 1 / 2) * g(d)) * pressure


def d_dot_ld2(d, temperature, pressure):
    return 15 * fun_aux.dv(temperature) * material.OMEGA / (
            constants.BOLTZMANN * temperature * np.power(material.R,
                                                         2)) * np.power(d,
                                                                        2) * np.power(
        (1 - d) / (5 * d), 1 / 3) / (
                   1 - 5 * np.power((1 - d) / (5 * d), 2 / 3)) * pressure


def d_dot_nhc1(d, temperature, pressure):
    if material.G < 2 * fun_aux.x(d):
        return 4.0 * material.OMEGA / (
                constants.BOLTZMANN * temperature * np.power(material.G,
                                                             2)) * np.power(
            constants.D0 / d, 1 / 3) * (
                       fun_aux.dv(temperature) + np.pi * fun_aux.db(
                   temperature) / material.G) * pressure
    else:
        return 0


def d_dot_nhc2(d, temperature, pressure):
    if material.G < 2 * fun_aux.x(d):
        return 27 * material.OMEGA / (
                constants.BOLTZMANN * temperature * np.power(material.G,
                                                             2)) * (
                       1 - d) * (
                       fun_aux.dv(temperature) + np.pi * fun_aux.db(
                   temperature) / material.G) * pressure
    else:
        return 0


def d_total_1(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    elif np.isnan(d):
        d = 1.0
    elif d <= constants.D0:
        d = constants.D0 + 1E-15
    return d_dot_plc1(d, temperature, pressure) + \
           d_dot_gbd1(d, temperature, pressure) + \
           d_dot_ld1(d, temperature, pressure) + \
           d_dot_nhc1(d, temperature, pressure)


def d_total_2(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    elif np.isnan(d):
        d = 1.0
    elif d <= constants.D0:
        d = constants.D0 + 1E-15
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
    return d_dot_plc2(d, temperature, pressure) + \
           d_dot_gbd2(d, temperature, pressure) + \
           d_dot_ld2(d, temperature, pressure) + \
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

    # if d == 0.9:
    #     d = 0.9+1.0E-3

    return 0.5 * (1 + (np.exp(75 * (d - constants.DBREAK)) - 1) / (
            np.exp(75 * (d - constants.DBREAK)) + 1))


def d_total(t, d, temperature, pressure):
    return (1 - smoothing(d)) * d_total_1(d, temperature,
                                          pressure) + \
           smoothing(d) * d_total_2(d, temperature, pressure)

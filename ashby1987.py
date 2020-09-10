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


def c(d):
    return (1 - constants.D0) / (d - constants.D0)


def d_crp(temperature):
    return 1E-6 * np.exp(
        -material.QCRP / (constants.IDEALGASCONSTANT * temperature)) * (
                   material.TM / temperature - 2)


def pressure_i(d):
    if d >= constants.DC:
        return constants.P0
    else:
        return (1 - constants.DC) / (1 - d) * d / constants.DC * constants.P0


def d_yield1(temperature, pressure):
    return np.power(
        (1 - constants.D0) / (1.3 * material.SIGMAY0) * pressure + np.power(
            constants.D0, 3), 1 / 3)


def d_yield2(temperature, pressure):
    return 1 - np.exp(-3 / 2 * pressure / material.SIGMAY0)


def d_dot_plc1(d, temperature, pressure):
    return 3.1 / np.power(c(d), 1 / 2) * d_crp(temperature) * d * np.power(
        c(d) * (pressure - constants.P0) / (
                3 * material.SIGMAY0 * np.power(d, 2)), material.n)


def d_dot_plc2(d, temperature, pressure):
    return 1.5 * d_crp(temperature) * d * (1 - d) * np.power(
        1.5 / material.n * (pressure - pressure_i(d)) / material.SIGMAY0 / (
                1 - np.power(1 - d, 1 / material.n)), material.n)


def d_dot_bd1(d, temperature, pressure):
    return 43 * np.power(c(d), 2) * fun_aux.db(temperature) / np.power(
        material.R, 3) * (pressure - constants.P0) * material.OMEGA / (
                   constants.BOLTZMANN * temperature)


def d_dot_bd2(d, temperature, pressure):
    return 270 * np.power(1 - d, 1 / 2) * fun_aux.db(temperature) / np.power(
        material.R, 3) * (pressure - pressure_i(d)) * material.OMEGA / (
                   constants.BOLTZMANN * temperature)


def d_dot_vd1(d, temperature, pressure):
    return 43 * c(d) * (1 - constants.D0) * fun_aux.dv(temperature) / np.power(
        material.R, 2) * (pressure - pressure_i(d)) * material.OMEGA / (
                   constants.BOLTZMANN * temperature)


def d_dot_vd2(d, temperature, pressure):
    return 270 * np.power(1 - d, 1 / 2) * np.power((1 - d) / 6,
                                                   1 / 3) * fun_aux.dv(
        temperature) / np.power(material.R, 2) * (
                   pressure - pressure_i(d)) * material.OMEGA / (
                   constants.BOLTZMANN * temperature)


def d_total_without_nhc1(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
    return d_dot_plc1(d, temperature, pressure) + \
           d_dot_bd1(d, temperature, pressure) + \
           d_dot_vd1(d, temperature, pressure)


def d_total_without_nhc2(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
    return d_dot_plc2(d, temperature, pressure) + \
           d_dot_bd2(d, temperature, pressure) + \
           d_dot_vd2(d, temperature, pressure)


def d_total_with_nhc1(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
    return d_dot_plc1(d, temperature, pressure) + \
           d_dot_bd1(d, temperature, pressure) + \
           d_dot_vd1(d, temperature, pressure)


def d_total_with_nhc2(d, temperature, pressure):
    if d > 1.0:
        d = 1.0
    # if d < constants.D0:
    #     d = constants.D0+1.0E-6
    return d_dot_plc2(d, temperature, pressure) + \
           d_dot_bd2(d, temperature, pressure) + \
           d_dot_vd2(d, temperature, pressure)


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


def d_total_without_nhc(t, d, temperature, pressure):
    return (1 - smoothing(d)) * d_total_without_nhc1(d, temperature,
                                                     pressure) + \
           smoothing(d) * d_total_without_nhc2(d, temperature, pressure)


def d_total_with_nhc(t, d, temperature, pressure):
    return (1 - smoothing(d)) * d_total_with_nhc1(d, temperature, pressure) + \
           smoothing(d) * d_total_with_nhc2(d, temperature, pressure)

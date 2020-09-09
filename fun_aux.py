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
import material


def mu(temperature):
    return material.MU0 * (
                1 + (temperature - 300) / material.TM * material.BETA)


def sigmay(temperature):
    return material.SIGMAY0 * (
                1 + (temperature - 300) / material.TM * material.ALPHA)


def dv(temperature):
    return material.D0V * np.exp(
        -material.QV / (constants.IDEALGASCONSTANT * temperature))


def db(temperature):
    return material.DD0B * np.exp(
        -material.QB / (constants.IDEALGASCONSTANT * temperature))


def epssig(temperature):
    return material.A * material.B * dv(temperature) / (
            constants.BOLTZMANN * temperature * np.power(mu(temperature),
                                                         material.n - 1))


def x(d):
    if d < constants.D0:
        d = constants.D0
    if d > 1:
        d = 1
    return np.power(3, -1 / 2) * np.power(
        (d - constants.D0) / (1 - constants.D0), 1 / 2) * material.R


def p_eff(d, pressure):
    return pressure * (1 - constants.D0) / (
                np.power(d, 2) * (d - constants.D0))


def rho(d):
    return material.R * (d - constants.D0)


def r(d):
    if d > 1:
        d = 1
    return material.R * np.power((1 - d) / 6, 1 / 3)


def ramps():
    T = np.array([6, 2])
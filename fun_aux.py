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
        -material.QV / (constants.GASCONSTANT * temperature))


def db(temperature):
    return material.DD0B * np.exp(
        -material.QB / (constants.GASCONSTANT * temperature))


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

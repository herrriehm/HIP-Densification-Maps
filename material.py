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

R = 50.0E-6  # m
G = 100.0E-6  # m


# # Copper
# OMEGA = 1.18E-29
# B = 2.56E-10
# TM = 1356
# SIGMAY = 50.0E6
# SIGMAY0 = 0
# ALPHA = 0
# MU0 = 4.21E10
# BETA = -0.54
# D0V = 2.0E-5
# QV = 197.0E3
# DD0B = 5.0E-15
# QB = 104.0E3
# n = 4.8
# A = 7.4E5
#
# R = 75.0E-6
# G = 10.0E-6
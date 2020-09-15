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
SIGMAY = 200.0E6  # Pa; from Helle 1985
SIGMAY0 = 200.0E6  # Pa; from Redouani 2015
ALPHA = -0.85  # from Redouani 2015
# SIGMAY0 = 9.1E8  # Pa; from Arzt 1983
# ALPHA = -4.5  # from Arzt 1983
MU0 = 8.1E10  # Pa
BETA = -0.85
D0V = 3.7E-5  # m^2 s^-1
QV = 280.0E3  # J mol^-1
QCRP = 280.0E3  # J mol^-1
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
# SIGMAY0 = 55.0E6  # from Redouani 2015
# ALPHA = -0.54  # from Redouani 2015
# MU0 = 4.21E10
# BETA = -0.54
# D0V = 2.0E-5
# QV = 197.0E3
# # DD0B = 5.0E-15
# DD0B = 5.12E-15  # from Redouani 2015
# QB = 104.0E3
# n = 4.8
# A = 7.4E5
#
# R = 100.0E-6
# G = 10.0E-6

# # Tungsten
# OMEGA = 1.59E-29
# B = 2.74E-10
# TM = 3680
# SIGMAY0 = 197.0E6
# ALPHA = -0.38
# MU0 = 3.45E11
# BETA = -0.38
# D0V = 5.6E-4
# QV = 585.0E3
# DD0B = 5.48E-13
# QB = 378.0E3
# n = 4.7
# A = 2.7E9
#
# R = 5.3E-6
# G = 10.0E-6

# # Tantalum
# OMEGA = 1.8E-29
# B = 2.86E-10
# TM = 3271
# SIGMAY0 = 600.5E6
# ALPHA = -0.23810
# MU0 = 6.12E10
# BETA = -0.42
# D0V = 1.2E-5
# QV = 413.0E3
# DD0B = 5.7E-14
# QB = 280.0E3
# n = 4.2
# A = 7.5E5
#
# R = 5.242E-6
# G = 1.0484E-05

# # 316L
# OMEGA = 1.21E-29
# B = 2.58E-10
# TM = 1680
# SIGMAY0 = 550.0E6  # 363.6E6  # fitted with Matlab for CvN Figure 5-5
# ALPHA = -0.85  # taken from Tool Steel data
# MU0 = 8.1E10
# BETA = -0.85
# D0V = 3.7E-5
# QV = 280.0E3
# QCRP = 280.0E3
# DD0B = 2.0E-13
# QB = 167.0E3
# n = 7.9
# A = 1.0E10
#
# R = 143.0E6  # m # taken with WPD from CvN Figure 6-1
# G = 10.0E-6  # m

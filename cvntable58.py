import numpy as np
# No 1
# TArray = np.array(
#     [25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25]) + 273
#
# PArray = np.array(
#     [0.1, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]) * 1000000

# # No 2
# TArray = np.array(
#     [25, 50, 75, 100, 125, 150, 175, 200, 200, 200, 200, 200, 200, 200, 200,
#      200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 217, 235,
#      252, 270, 287, 300, 300, 300, 300, 300]) + 273
#
# PArray = np.array(
#     [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 150, 150, 150, 150, 150]) * 1000000

# # No 3
# TArray = np.array(
#     [25, 50, 75, 100, 125, 150, 175, 200, 200, 200, 200, 200, 200, 200, 200,
#      200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 217, 235,
#      252, 270, 287, 305, 322, 340, 357, 375, 392, 410, 427, 445, 462, 480, 497,
#      515, 532, 550, 567, 585, 600, 600, 600, 600, 600]) + 273
#
# PArray = np.array(
#     [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 150, 150, 150, 150, 150]) * 1000000

# # No 4
# TArray = np.array(
#     [25, 50, 75, 100, 125, 150, 175, 200, 200, 200, 200, 200, 200, 200, 200,
#      200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 217, 235,
#      252, 270, 287, 305, 322, 340, 357, 375, 392, 410, 427, 445, 462, 480, 497,
#      515, 532, 550, 567, 585, 602, 620, 637, 655, 672, 690, 707, 724, 742, 759,
#      777, 794, 812, 829, 847, 864, 882, 900, 900, 900, 900, 900]) + 273
#
# PArray = np.array(
#     [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
#      0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 150, 150, 150, 150, 150]) * 1000000

# # No 5
TArray = np.array(
    [25, 50, 75, 100, 125, 150, 175, 200, 200, 200, 200, 200, 200, 200, 200,
     200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 217, 235,
     252, 270, 287, 305, 322, 340, 357, 375, 392, 410, 427, 445, 462, 480, 497,
     515, 532, 550, 567, 585, 602, 620, 637, 655, 672, 690, 707, 724, 742, 759,
     777, 794, 812, 829, 847, 864, 882, 899, 917, 934, 952, 969, 987, 1004,
     1022, 1039, 1057, 1074, 1092, 1109, 1125, 1125, 1125, 1125, 1125]) + 273

PArray = np.array(
    [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
     0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
     0.1, 0.1, 0.1, 0.1, 0.1, 30, 30, 30, 30, 30]) * 1000000

# HIP process parameters from Reardon 1998 "HIP 9.0 User's Guide"
# this should really be encapsulated in some nice small functions

import numpy as np

TArray = np.array(
    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38,
     39, 41, 42, 44, 45, 47, 48, 50, 51, 53, 54, 56, 57, 59, 60, 62, 63, 65,
     65, 66, 66, 67, 67, 68, 68, 68, 69, 69, 70, 70, 71, 71, 71, 72, 72, 73,
     73, 74, 74, 75, 75, 75, 76, 80, 89, 98, 107, 116, 125, 134, 143, 152, 161,
     170, 179, 188, 197, 206, 215, 224, 233, 242, 251, 260, 269, 278, 287, 296,
     305, 314, 323, 332, 341, 349, 358, 367, 376, 385, 394, 403, 412, 421, 430,
     439, 448, 457, 466, 475, 484, 493, 502, 511, 520, 529, 538, 547, 556, 565,
     574, 583, 592, 601, 610, 619, 628, 637, 646, 655, 664, 673, 682, 691, 700,
     709, 718, 727, 735, 744, 753, 762, 771, 780, 789, 798, 807, 816, 825, 834,
     843, 852, 861, 870, 879, 888, 897, 906, 915, 924, 933, 942, 951, 960, 969,
     978, 987, 996, 1005, 1014, 1023, 1032, 1041, 1050, 1059, 1068, 1077, 1086,
     1095, 1104, 1110, 1117, 1123, 1130, 1137, 1143, 1150, 1150, 1150, 1150,
     1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150,
     1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150,
     1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150,
     1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150,
     1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1150, 1105, 1060, 1014,
     969, 924, 879, 834, 788, 743, 698, 653, 608, 562, 517, 472, 427, 382, 336,
     291, 246, 201, 156, 110, 65, 20]) + 273

PArray = np.array(
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 14, 28, 41, 55, 68, 74, 79, 85, 90,
     96, 101, 107, 112, 118, 123, 129, 134, 140, 145, 151, 156, 162, 167, 173,
     178, 184, 189, 195, 202, 209, 217, 224, 231, 238, 245, 253, 260, 267, 274,
     282, 289, 296, 303, 310, 318, 325, 332, 339, 346, 354, 361, 368, 375, 383,
     390, 397, 404, 411, 419, 426, 433, 440, 447, 455, 462, 469, 476, 484, 491,
     498, 505, 512, 520, 527, 534, 539, 545, 550, 555, 560, 566, 571, 576, 582,
     587, 592, 598, 603, 608, 613, 619, 624, 629, 635, 640, 645, 650, 656, 661,
     666, 672, 677, 682, 688, 693, 698, 703, 709, 714, 719, 725, 730, 735, 741,
     746, 751, 756, 762, 767, 772, 778, 783, 788, 793, 799, 804, 809, 815, 820,
     825, 831, 836, 841, 846, 852, 857, 862, 868, 873, 878, 883, 889, 894, 899,
     905, 910, 915, 921, 926, 931, 936, 942, 947, 952, 958, 963, 968, 974, 979,
     984, 989, 995, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
     1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 960, 920, 880,
     840, 800, 760, 720, 680, 640, 600, 560, 520, 481, 441, 401, 361, 321, 281,
     241, 201, 161, 121, 81, 41, 1]) * 100000
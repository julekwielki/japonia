import math

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from Av_lin_fits_class import Fitted_one_log
from Av_lin_fits_class import Fitted


time00 = [33, 38, 47, 48, 55, 68, 70, 73, 91]
survival00 = [0.962, 0.923, 0.885, 0.846, 0.806, 0.761, 0.714, 0.618, 0.515]
std00 = [0.0377, 0.0523, 0.0627, 0.0708, 0.078, 0.0856, 0.0925, 0.1017, 0.1266]

time01 = [85, 86, 88, 96]
survival01 = [0.952, 0.857, 0.807, 0.733]
std01 = [0.0465, 0.0764, 0.0869, 0.1055]

time30 = [30, 32, 33, 36, 37, 43, 45, 46, 51, 53, 55, 60]
survival30 = [0.96, 0.92, 0.88, 0.8, 0.76, 0.715, 0.671, 0.626, 0.578, 0.525, 0.473, 0.414]
std30 = [0.0392, 0.0543, 0.065, 0.08, 0.0854, 0.0913, 0.096, 0.0994, 0.1028, 0.106, 0.1076, 0.1092]

time31 = [35, 36, 52, 64, 65, 66, 68, 69, 72, 76, 77, 80, 82]
survival31 = [0.967, 0.933, 0.896, 0.815, 0.772, 0.729, 0.686, 0.643, 0.594, 0.534, 0.475, 0.396, 0.297]
std31 = [0.0328, 0.0455, 0.057, 0.0755, 0.0828, 0.0886, 0.0932, 0.0967, 0.1012, 0.107, 0.1104, 0.117, 0.1226]

time70 = [10, 13, 16, 26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80]
survival70 = [0.962, 0.923, 0.885, 0.808, 0.769, 0.731, 0.692, 0.654, 0.613, 0.572, 0.531, 0.49, 0.45, 0.409, 0.368,
              0.315, 0.263, 0.131]
std70 = [0.0377, 0.0523, 0.0627, 0.0773, 0.0826, 0.087, 0.0905, 0.0933, 0.096, 0.0979, 0.0991, 0.0995, 0.0993, 0.0983,
         0.0966, 0.096, 0.0933, 0.1039]

time71 = [12, 13, 30, 37, 39, 40, 42, 44, 45, 48, 55, 64, 81, 82]
survival71 = [0.964, 0.821, 0.784, 0.709, 0.672, 0.635, 0.597, 0.56, 0.523, 0.448, 0.392, 0.336, 0.252, 0.126]
std71 = [0.0351, 0.0724, 0.0781, 0.0867, 0.0898, 0.0923, 0.0941, 0.0953, 0.096, 0.0957, 0.0988, 0.0993, 0.1041, 0.1032]


# """
def func(x, a, n):
    return np.log(a) + n * np.log(x)

t = [time00, time01, time30, time31, time70, time71]
surv = [survival00, survival01, survival30, survival31, survival70, survival71]
stdev = [std00, std01, std30, std31, std70, std71]
name = ["t0 p0", "t0 p1", "t3 p0", "t3 p1", "t7 p0", "t7 p1"]
col = ['blue', 'orange', 'green', 'red', 'plum', 'brown']


text = ""
for i in range(len(t)):
    xx = t[i]
    xx2 = [x*7/365 for x in xx]
    sd = stdev[i]
    data = [1-x for x in surv[i]]
    yy = [np.log(-np.log(1 - y)) for y in data]
    sd_log = [np.abs(sd[i] / (np.log(1 - data[i]) * (1 - data[i]))) for i in range(len(sd))]
    A = Fitted_one_log(xx2, data, sd)
    popt, pcov, perr, r_squared = A.fit_both_sd()

    plt.errorbar(popt[0], popt[1], xerr=perr[0], yerr=perr[1], fmt='o', ecolor=col[i], color=col[i], label=name[i])

N_non_mut_x_weeks = [21, 38, 59, 76]
N_non_mut_x_years = [x * 7 / 365 for x in N_non_mut_x_weeks]
N_non_mut_sd = [0.0377, 0.0545, 0.0685, 0.098]
N_non_mut_no_canc = [0.962, 0.92, 0.874, 0.794]
N_non_mut_canc = [1 - x for x in N_non_mut_no_canc]

A = Fitted_one_log(N_non_mut_x_years, N_non_mut_canc, N_non_mut_sd)
popt1, pcov1, perr1, r_squared1 = A.fit_both_sd()
plt.errorbar(popt1[0], popt1[1], xerr=perr1[0], yerr=perr1[1], fmt='o', label="BRCA1 normal")

N_mut_x_weeks = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
N_mut_x_years = [x * 7 / 365 for x in N_mut_x_weeks]
N_mut_sd = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
N_mut_no_canc = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
N_mut_canc = [1 - x for x in N_mut_no_canc]

A = Fitted_one_log(N_mut_x_years, N_mut_canc, N_mut_sd)
popt2, pcov2, perr2, r_squared2 = A.fit_both_sd()
plt.errorbar(popt2[0], popt2[1], xerr=perr2[0], yerr=perr2[1], fmt='o', ecolor='indigo', color='indigo', label="BRCA1 mutation")

plt.legend()
plt.xlabel("parameter a")
plt.ylabel("parameter k")
# plt.xscale("log")
# plt.ylabel("log(-log(1 - x))")
plt.show()


time00 = [33, 38, 47, 48, 55, 68, 70, 73, 91]
survival00 = [0.962, 0.923, 0.885, 0.846, 0.806, 0.761, 0.714, 0.618, 0.515]
std00 = [0.0377, 0.0523, 0.0627, 0.0708, 0.078, 0.0856, 0.0925, 0.1017, 0.1266]

time01 = [85, 86, 88, 96]
survival01 = [0.952, 0.857, 0.807, 0.733]
std01 = [0.0465, 0.0764, 0.0869, 0.1055]

time30 = [30, 32, 33, 36, 37, 43, 45, 46, 51, 53, 55, 60]
survival30 = [0.96, 0.92, 0.88, 0.8, 0.76, 0.715, 0.671, 0.626, 0.578, 0.525, 0.473, 0.414]
std30 = [0.0392, 0.0543, 0.065, 0.08, 0.0854, 0.0913, 0.096, 0.0994, 0.1028, 0.106, 0.1076, 0.1092]

time31 = [35, 36, 52, 64, 65, 66, 68, 69, 72, 76, 77, 80, 82]
survival31 = [0.967, 0.933, 0.896, 0.815, 0.772, 0.729, 0.686, 0.643, 0.594, 0.534, 0.475, 0.396, 0.297]
std31 = [0.0328, 0.0455, 0.057, 0.0755, 0.0828, 0.0886, 0.0932, 0.0967, 0.1012, 0.107, 0.1104, 0.117, 0.1226]

T_20w_time70 = [26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80]
T_20w_survival70 = [0.913, 0.87, 0.826, 0.783, 0.739, 0.693, 0.647, 0.601, 0.554, 0.508, 0.462, 0.416, 0.356, 0.297,
                    0.148]
T_20w_std70 = [0.0588, 0.0702, 0.079, 0.086, 0.0916, 0.0968, 0.1008, 0.1036, 0.1054, 0.1063, 0.1062, 0.1051, 0.1056,
               0.1033, 0.117]

T_20w_time71 = [30, 37, 39, 40, 42, 44, 45, 48, 55, 64, 81, 82]
T_20w_survival71 = [0.955, 0.864, 0.818, 0.773, 0.727, 0.682, 0.636, 0.545, 0.477, 0.409, 0.307, 0.153]
T_20w_std71 = [0.0444, 0.0732, 0.0822, 0.0893, 0.095, 0.0993, 0.1026, 0.1062, 0.1127, 0.1154, 0.1238, 0.1249]

# """
def func(x, a, n):
    return np.log(a) + n * np.log(x)

t = [time00, time01, time30, time31, T_20w_time70, T_20w_time71]
surv = [survival00, survival01, survival30, survival31, T_20w_survival70, T_20w_survival71]
stdev = [std00, std01, std30, std31, T_20w_std70, T_20w_std71]
name = ["t0 p0", "t0 p1", "t3 p0", "t3 p1", "t7 p0 late", "t7 p1 late"]
col = ['blue', 'orange', 'green', 'red', 'plum', 'brown']


text = ""
for i in range(len(t)):
    xx = t[i]
    xx2 = [x*7/365 for x in xx]
    sd = stdev[i]
    data = [1-x for x in surv[i]]
    yy = [np.log(-np.log(1 - y)) for y in data]
    sd_log = [np.abs(sd[i] / (np.log(1 - data[i]) * (1 - data[i]))) for i in range(len(sd))]
    A = Fitted_one_log(xx2, data, sd)
    popt, pcov, perr, r_squared = A.fit_both_sd()

    plt.errorbar(popt[0], popt[1], xerr=perr[0], yerr=perr[1], fmt='o', ecolor=col[i], color=col[i], label=name[i])

N_non_mut_x_weeks = [21, 38, 59, 76]
N_non_mut_x_years = [x * 7 / 365 for x in N_non_mut_x_weeks]
N_non_mut_sd = [0.0377, 0.0545, 0.0685, 0.098]
N_non_mut_no_canc = [0.962, 0.92, 0.874, 0.794]
N_non_mut_canc = [1 - x for x in N_non_mut_no_canc]

A = Fitted_one_log(N_non_mut_x_years, N_non_mut_canc, N_non_mut_sd)
popt1, pcov1, perr1, r_squared1 = A.fit_both_sd()
plt.errorbar(popt1[0], popt1[1], xerr=perr1[0], yerr=perr1[1], fmt='o', label="BRCA1 normal")

N_mut_x_weeks = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
N_mut_x_years = [x * 7 / 365 for x in N_mut_x_weeks]
N_mut_sd = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
N_mut_no_canc = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
N_mut_canc = [1 - x for x in N_mut_no_canc]

A = Fitted_one_log(N_mut_x_years, N_mut_canc, N_mut_sd)
popt2, pcov2, perr2, r_squared2 = A.fit_both_sd()
plt.errorbar(popt2[0], popt2[1], xerr=perr2[0], yerr=perr2[1], fmt='o', ecolor='indigo', color='indigo', label="BRCA1 mutation")

plt.legend()
plt.xlabel("parameter a")
plt.ylabel("parameter k")
# plt.xscale("log")
# plt.ylabel("log(-log(1 - x))")
plt.show()


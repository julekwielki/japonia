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

time70 = [26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80]
survival70 = [0.913, 0.87, 0.826, 0.783, 0.739, 0.693, 0.647, 0.601, 0.554, 0.508, 0.462, 0.416, 0.356, 0.297, 0.148]
std70 = [0.0588, 0.0702, 0.079, 0.086, 0.0916, 0.0968, 0.1008, 0.1036, 0.1054, 0.1063, 0.1062, 0.1051, 0.1056, 0.1033, 0.117]

time71 = [30, 37, 39, 40, 42, 44, 45, 48, 55, 64, 81, 82]
survival71 = [0.955, 0.864, 0.818, 0.773, 0.727, 0.682, 0.636, 0.545, 0.477, 0.409, 0.307, 0.153]
std71 = [0.0444, 0.0732, 0.0822, 0.0893, 0.095, 0.0993, 0.1026, 0.1062, 0.1127, 0.1154, 0.1238, 0.1249]

t = [time00, time01, time30, time31, time70, time71]
surv = [survival00, survival01, survival30, survival31, survival70, survival71]
stdev = [std00, std01, std30, std31, std70, std71]
name = ["t0 p0", "t0 p1", "t3 p0", "t3 p1", "t7 p0 late", "t7 p1 late"]
col = ['blue', 'red', 'orange', 'plum', 'green', 'brown']

text = ""

fig, ax = plt.subplots(2, 1, figsize=(8, 15))

for i in range(len(t)):
    xx = t[i]
    xx2 = [x*7/365 for x in xx]
    sd = stdev[i]
    data = [1-x for x in surv[i]]

    A = Fitted(xx2, data, sd)
    A.fit_both_sd()

    dane1, dane2 = A.toasty()

    ax[0].errorbar(xx2, data, sd, fmt='o',  ecolor=col[i], color=col[i], label=name[i])

    ax[0].plot(xx2, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in xx2], '--', color=col[i])

    ax[1].errorbar(xx2, [np.log(-np.log(1 - y)) for y in data],
                   [np.abs(sd[i] / (np.log(1 - data[i]) * (1 - data[i]))) for i in range(len(sd))], fmt='o',  ecolor=col[i], color=col[i], label=name[i])

    ax[1].plot(xx2, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in xx2], '--', color=col[i])

    # print(name[i])
    # print(popt)
    # print(perr)
    # print("r^2", r_squared)

    # text = text + name[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{r_squared:.3f}' + "\n"



ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")
ax[0].legend(loc="upper left")
ax[1].legend(loc="upper left")
plt.show()
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from Av_lin_fits_class import Fitted_one_log

x_contr = [43, 46, 50, 52, 61, 63, 64, 65, 67, 69, 70, 71, 73, 74, 75, 77, 78, 79, 82, 83, 84, 85, 89]
y_contr = [0.987, 0.974, 0.961, 0.934, 0.92, 0.907, 0.893, 0.878, 0.864, 0.85, 0.835, 0.821, 0.789, 0.757, 0.74, 0.706,
           0.689, 0.671, 0.653, 0.633, 0.614, 0.554, 0.532]
sd_contr = [0.0131, 0.0184, 0.0223, 0.0284, 0.0312, 0.0337, 0.0359, 0.038, 0.04, 0.0418, 0.0435, 0.0451, 0.0486,
              0.0518, 0.0532, 0.0559, 0.0572, 0.0584, 0.0597, 0.0609, 0.0622, 0.0649, 0.066]

x_3 = [21, 42, 47, 53, 79, 84]
y_3 = [0.958, 0.915, 0.818, 0.77, 0.693, 0.607]
sd_3 = [0.0408, 0.0577, 0.0825, 0.0906, 0.1095, 0.1255]

x_6 = [43, 50, 58, 72, 74, 89]
y_6 = [0.958, 0.917, 0.873, 0.818, 0.764, 0.688]
sd_6 = [0.0408, 0.0564, 0.0686, 0.0832, 0.0939, 0.1113]

x_12 = [37, 43, 52, 53, 58, 63, 66, 74, 75, 84]
y_12 = [0.958, 0.875, 0.833, 0.789, 0.746, 0.702, 0.655, 0.605, 0.554, 0.462]
sd_12 = [0.0408, 0.0675, 0.0761, 0.0838, 0.0899, 0.0947, 0.0993, 0.1036, 0.1065, 0.1224]

x_24 = [18, 22, 33, 36, 51, 52, 53, 58]
y_24 = [0.958, 0.917, 0.875, 0.833, 0.789, 0.743, 0.693, 0.644]
sd_24 = [0.0408, 0.0564, 0.0675, 0.0761, 0.0838, 0.0908, 0.0973, 0.1022]

x_60 = [12, 29, 32, 35, 37, 38, 44, 45, 50, 55]
y_60 = [0.957, 0.913, 0.87, 0.783, 0.652, 0.609, 0.565, 0.522, 0.435, 0.386]  # data
sd_60 = [0.0425, 0.0588, 0.0702, 0.086, 0.0993, 0.1018, 0.1034, 0.1042, 0.1034, 0.1026]

XX = [x_contr, x_3, x_6, x_12, x_24, x_60]
YY = [y_contr, y_3, y_6, y_12, y_24, y_60]
SD = [sd_contr, sd_3, sd_6, sd_12, sd_24, sd_60]

dose_rates = [0, 3, 6, 12, 24, 60]

def function(x, a, b, c):
    return np.exp(- a * x**2 * np.exp(-b * x)) + c * x


def function_av(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))

data = [1-x for x in y_contr]
xx = x_contr
sd = sd_contr
A = Fitted_one_log(x_contr, data, sd_contr)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2']]

fig, ax = plt.subplots(2, 1, figsize=(8, 15))
ax[0].errorbar(xx, data, yerr=sd, fmt='o')
ax[0].plot(xx, [function_av(x, aaa[0], aaa[2]) for x in xx], '--')

# ax[1].scatter(xx, A.y, label=name0[i], color=col0[i])
ax[1].errorbar(xx, A.y, A.sd, fmt='o')
ax[1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x], '--')
ax[1].set_xscale('log')
plt.show()
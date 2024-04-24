import Avrami_data
from Av_lin_fits_class import Fitted
import matplotlib.pyplot as plt
import numpy as np

"""T_time70 = [10, 13, 16, 26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80, 12, 13, 30, 37, 39, 40, 42, 44, 45, 48, 55, 64, 81, 82]
T_survival70 = [0.962, 0.923, 0.885, 0.808, 0.769, 0.731, 0.692, 0.654, 0.613, 0.572, 0.531, 0.49, 0.45, 0.409, 0.368,
                0.315, 0.263, 0.131, 0.964, 0.821, 0.784, 0.709, 0.672, 0.635, 0.597, 0.56, 0.523, 0.448, 0.392, 0.336, 0.252, 0.126]
T_std70 = [0.0377, 0.0523, 0.0627, 0.0773, 0.0826, 0.087, 0.0905, 0.0933, 0.096, 0.0979, 0.0991, 0.0995, 0.0993, 0.0983,
           0.0966, 0.096, 0.0933, 0.1039, 0.0351, 0.0724, 0.0781, 0.0867, 0.0898, 0.0923, 0.0941, 0.0953, 0.096, 0.0957, 0.0988, 0.0993, 0.1041,
           0.1032]

x1 = [x * 7/365 for x in T_time70]
y1 = [1-x for x in T_survival70]
sd1 = T_std70

A = Fitted(x1, y1, sd1)
A.fit_both_sd()

dane1, dane2 = A.toasty()

T_time70 = [10, 13, 16, 26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80]
T_survival70 = [0.962, 0.923, 0.885, 0.808, 0.769, 0.731, 0.692, 0.654, 0.613, 0.572, 0.531, 0.49, 0.45, 0.409, 0.368,
                0.315, 0.263, 0.131]
T_std70 = [0.0377, 0.0523, 0.0627, 0.0773, 0.0826, 0.087, 0.0905, 0.0933, 0.096, 0.0979, 0.0991, 0.0995, 0.0993, 0.0983,
           0.0966, 0.096, 0.0933, 0.1039]
x2 = [x * 7/365 for x in T_time70]
T_time71 = [12, 13, 30, 37, 39, 40, 42, 44, 45, 48, 55, 64, 81, 82]
T_survival71 = [0.964, 0.821, 0.784, 0.709, 0.672, 0.635, 0.597, 0.56, 0.523, 0.448, 0.392, 0.336, 0.252, 0.126]
T_std71 = [0.0351, 0.0724, 0.0781, 0.0867, 0.0898, 0.0923, 0.0941, 0.0953, 0.096, 0.0957, 0.0988, 0.0993, 0.1041,
           0.1032]
print(dane2)
fig, ax = plt.subplots(2, 1, figsize=(8, 15))
ax[0].scatter(x1, y1)
ax[0].errorbar(x1, y1, sd1, fmt='o')  

ax[0].plot(x2, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--') 

ax[1].scatter(x1, [np.log(-np.log(1 - y)) for y in y1])
ax[1].errorbar(x1, [np.log(-np.log(1 - y)) for y in y1],[np.abs(sd1[i] / (np.log(1 - y1[i]) * (1 - y1[i]))) for i in range(len(sd1))], fmt='o') 

ax[1].plot(x2, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--')
ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")
plt.show() """

# po 20

T_time70 = [26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80, 30, 37, 39, 40, 42, 44, 45, 48, 55, 64, 81, 82]
T_survival70 = [0.913, 0.87, 0.826, 0.783, 0.739, 0.693, 0.647, 0.601, 0.554, 0.508, 0.462, 0.416, 0.356, 0.297,
                    0.148, 0.955, 0.864, 0.818, 0.773, 0.727, 0.682, 0.636, 0.545, 0.477, 0.409, 0.307, 0.153]
T_std70 = [0.0588, 0.0702, 0.079, 0.086, 0.0916, 0.0968, 0.1008, 0.1036, 0.1054, 0.1063, 0.1062, 0.1051, 0.1056,
               0.1033, 0.117, 0.0444, 0.0732, 0.0822, 0.0893, 0.095, 0.0993, 0.1026, 0.1062, 0.1127, 0.1154, 0.1238, 0.1249]


x1 = [x * 7/365 for x in T_time70]
y1 = [1-x for x in T_survival70]
sd1 = T_std70

A = Fitted(x1, y1, sd1)
A.fit_both_sd()

dane1, dane2 = A.toasty()

T_time70 = [26, 30, 31, 32, 35, 39, 40, 42, 43, 47, 50, 63, 66, 67, 80]
x2 = [x * 7/365 for x in T_time70]

print(dane2)
fig, ax = plt.subplots(2, 1, figsize=(8, 15))
ax[0].scatter(x1, y1)
ax[0].errorbar(x1, y1, sd1, fmt='o')  # """

ax[0].plot(x2, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--') # """

ax[1].scatter(x1, [np.log(-np.log(1 - y)) for y in y1])
ax[1].errorbar(x1, [np.log(-np.log(1 - y)) for y in y1],[np.abs(sd1[i] / (np.log(1 - y1[i]) * (1 - y1[i]))) for i in range(len(sd1))], fmt='o')  # """

ax[1].plot(x2, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--')
ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")
plt.show()
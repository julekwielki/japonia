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

"""

def func(x, c, a, n):
    return c * (1 - np.exp(-a * np.power(x, n)))


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

    popt, pcov = curve_fit(func, xx2, data, sigma=sd, bounds=([0., 0, 0], [1, 5., 9.]))
    perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe

    residuals = data - func(xx2, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((data - np.mean(data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    plt.errorbar(xx, data, yerr=sd, fmt='o', ecolor=col[i], color=col[i])
    plt.scatter(xx, data, label=name[i], color=col[i])
    plt.plot(xx, [func(x, popt[0], popt[1], popt[2]) for x in xx2], '--', label='fit: n=%.1f' % (popt[2]), color=col[i])

    # print(name[i])
    # print(popt)
    # print(perr)
    # print("r^2", r_squared)

    # text = text + name[i] + "\nc(u(c))\ta(u(a))\tn(u(n))\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{popt[2]:.3f}' + " (" + f'{perr[2]:.3f}' + ")\t" + f'{r_squared:.3f}' + "\n"

    text = text + name[i] + "\nc\tu(c)\ta\tu(a)\tn\tu(n)\tr^2\n"
    text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{popt[2]:.3f}' + "\t" + f'{perr[2]:.3f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

plt.legend()
plt.xlabel("time (weeks)")
plt.ylabel("carcinoma probability")
plt.title("Parity fit - fit for c")

plt.show()
# """  # full fit with c

#"""

def func(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


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

    popt, pcov = curve_fit(func, xx2, data, sigma=sd, bounds=([0, 0], [1., 10.]))
    perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe

    residuals = data - func(xx2, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((data - np.mean(data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    plt.errorbar(xx, data, yerr=sd, fmt='o', ecolor=col[i], color=col[i])
    plt.scatter(xx, data, label=name[i], color=col[i])
    plt.plot(xx, [func(x, popt[0], popt[1]) for x in xx2], '--', label='fit: n=%.1f' % (popt[1]), color=col[i])

    # print(name[i])
    # print(popt)
    # print(perr)
    # print("r^2", r_squared)

    # text = text + name[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{r_squared:.3f}' + "\n"

    text = text + name[i] + "\na\tu(a)\tn\tu(n)\tr^2\n"
    text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

plt.legend()
plt.xlabel("time (weeks)")
plt.ylabel("carcinoma probability")
plt.title("Parity fit - c = 1")

plt.show()

# """  # fit for c = 1

"""

def func(x, c, a, n):
    return c * (1 - np.exp(-a * np.power(x, n)))


t = [time00, time31, time30, time71, time70]
surv = [survival00, survival31, survival30, survival71, survival70]
stdev = [std00, std31, std30, std71, std70]
name = ["t0 p0", "t3 p1", "t3 p0", "t7 p1", "t7 p0"]
col = ['blue', 'green', 'red', 'plum', 'brown']

xx = time01
xx2 = [x*7/365 for x in xx]
sd = std01
data = [1-x for x in survival01]

popt, pcov = curve_fit(func, xx2, data, sigma=sd, bounds=([0., 0, 0], [1, 1, 10.]))
perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe

residuals = data - func(xx2, *popt)
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((data - np.mean(data)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

plt.errorbar(xx, data, yerr=sd, fmt='o', ecolor='orange', color='orange')
plt.scatter(xx, data, label="t0 p1", color='orange')
plt.plot(xx, [func(x, popt[0], popt[1], popt[2]) for x in xx2], '--', label='fit: n=%.1f' % (popt[2]), color='orange')

# print("t0 p1")
# print(popt)
# print(perr)
# print("r^2", r_squared)

# text = "t0 p1\nc(u(c))\ta(u(a))\tn(u(n))\tr^2\n"
# text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{popt[2]:.3f}' + " (" + f'{perr[2]:.3f}' + ")\t" + f'{r_squared:.3f}' + "\n"

text = "t0 p1\nc\tu(c)\ta\tu(a)\tn\tu(n)\tr^2\n"
text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{popt[2]:.3f}' + "\t" + f'{perr[2]:.3f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

text = ""
for i in range(len(t)):
    xx = t[i]
    xx2 = [x*7/365 for x in xx]
    sd = stdev[i]
    data = [1-x for x in surv[i]]
    print(name[i])
    popt, pcov = curve_fit(func, xx2, data, sigma=sd, bounds=([0., 0, 0], [1, 5., 9.]))
    perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe

    residuals = data - func(xx2, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((data - np.mean(data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    plt.errorbar(xx, data, yerr=sd, fmt='o', ecolor=col[i], color=col[i])
    plt.scatter(xx, data, label=name[i], color=col[i])
    plt.plot(xx, [func(x, popt[0], popt[1], popt[2]) for x in xx2], '--', label='fit: n=%.1f' % (popt[2]), color=col[i])

    # print(name[i])
    # print(popt)
    # print(perr)
    # print("r^2", r_squared)

    # text = text + name[i] + "\nc(u(c))\ta(u(a))\tn(u(n))\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{popt[2]:.3f}' + " (" + f'{perr[2]:.3f}' + ")\t" + f'{r_squared:.3f}' + "\n"

    text = text + name[i] + "\nc\tu(c)\ta\tu(a)\tn\tu(n)\tr^2\n"
    text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{popt[2]:.3f}' + "\t" + f'{perr[2]:.3f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

plt.legend()
plt.xlabel("time (weeks)")
plt.ylabel("carcinoma probability")
plt.title("Parity fit - fit for c (p0t1)")

plt.show()
# """  # full fit with c, t0 p1 by hand

"""

def func(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


t = [time00, time31, time30, time71, time70]
surv = [survival00, survival31, survival30, survival71, survival70]
stdev = [std00, std31, std30, std71, std70]
name = ["t0 p0", "t3 p1", "t3 p0", "t7 p1", "t7 p0"]
col = ['blue', 'green', 'red', 'plum', 'brown']

xx = time01
xx2 = [x*7/365 for x in xx]
sd = std01
data = [1-x for x in survival01]

popt, pcov = curve_fit(func, xx2, data, sigma=sd, bounds=([0, 0], [1, 12.]))
perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe

residuals = data - func(xx2, *popt)
ss_res = np.sum(residuals ** 2)
ss_tot = np.sum((data - np.mean(data)) ** 2)
r_squared = 1 - (ss_res / ss_tot)

plt.errorbar(xx, data, yerr=sd, fmt='o', ecolor='orange', color='orange')
plt.scatter(xx, data, label="t0 p1", color='orange')
plt.plot(xx, [func(x, popt[0], popt[1]) for x in xx2], '--', label='fit: n=%.1f' % (popt[1]), color='orange')

# print("t0 p1")
# print(popt)
# print(perr)
# print("r^2", r_squared)

# text = "t0 p1\na(u(a))\tn(u(n))\tr^2\n"
# text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{r_squared:.3f}' + "\n"

text = "t0 p1\na\tu(a)\tn\tu(n)\tr^2\n"
text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

text = ""
for i in range(len(t)):
    xx = t[i]
    xx2 = [x*7/365 for x in xx]
    sd = stdev[i]
    data = [1-x for x in surv[i]]

    popt, pcov = curve_fit(func, xx2, data, sigma=sd, bounds=([0, 0], [1., 10.]))
    perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe

    residuals = data - func(xx2, *popt)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((data - np.mean(data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    plt.errorbar(xx, data, yerr=sd, fmt='o', ecolor=col[i], color=col[i])
    plt.scatter(xx, data, label=name[i], color=col[i])
    plt.plot(xx, [func(x, popt[0], popt[1]) for x in xx2], '--', label='fit: n=%.1f' % (popt[1]), color=col[i])

    # print(name[i])
    # print(popt)
    # print(perr)
    # print("r^2", r_squared)

    # text = text + name[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{r_squared:.3f}' + "\n"

    text = text + name[i] + "\na\tu(a)\tn\tu(n)\tr^2\n"
    text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

plt.legend()
plt.xlabel("time (weeks)")
plt.ylabel("carcinoma probability")
plt.title("Parity fit - c = 1 (t0 p1)")

plt.show()

# """  # fit for c = 1, t0 p1 by hand

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

    # plt.errorbar(xx2, yy, yerr=sd_log, fmt='o', ecolor=col[i], color=col[i])
    # plt.scatter(xx2, yy, label=name[i], color=col[i])
    # plt.plot(xx2, [func(x, popt[0], popt[1]) for x in xx2], '--', label='fit: n=%.1f' % (popt[1]), color=col[i])
    # print(name[i])
    # print(popt)
    # print(perr)
    # print("r^2", r_squared)

    # text = text + name[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{r_squared:.3f}' + "\n"

    # text = text + name[i] + "\na\tu(a)\tn\tu(n)\tr^2\n"
    # text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{r_squared:.3f}' + "\n"
    plt.errorbar(popt[0], popt[1], xerr=perr[0], yerr=perr[1], fmt='o', ecolor=col[i], color=col[i], label=name[i])

plt.legend()
# plt.xscale("log")
# plt.ylabel("log(-log(1 - x))")
plt.show()


# """  # log fit
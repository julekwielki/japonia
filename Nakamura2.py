import math

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

"""

def function(x, c, a, n):
    return c * (1 - np.exp(-a * np.power(x, n)))


x1 = [21, 38, 59, 76]
x1_scaled = [x*7/365 for x in x1]
sd1 = [0.0377, 0.0545, 0.0685, 0.098]
data = [0.962, 0.92, 0.874, 0.794]
data_can = [1-x for x in data]

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

popt, pcov = curve_fit(function, x1_scaled, data_can, sigma=sd1, bounds=([0, 0, 0], [1, 5., 6.]))
perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe
print(popt)
print(perr)

residuals = data_can - function(x1_scaled, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_can-np.mean(data_can))**2)
r_squared = 1 - (ss_res / ss_tot)
print("r^2", r_squared)
print(type(popt))

text = f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{popt[2]:.3f}' + " (" + f'{perr[2]:.3f}' + ")\t" + f'{r_squared:.3f}' + "\n"
print(text)

ax[0].plot(x1, [function(x, popt[0], popt[1], popt[2]) for x in x1_scaled], 'g--', label='fit: c=%.3f, a=%.5f, n=%.3f' % (popt[0], popt[1], popt[2]))
ax[0].errorbar(x1, data_can, yerr=sd1, fmt='o', label="data")
ax[0].legend()
ax[0].set_xlabel("time (weeks)")
ax[0].set_ylabel("carcinoma probability genotype 0 - BRCA1 +/+")


xx2 = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
xx22 = [x*7/365 for x in xx2]
data2 = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
sd2 = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
data_can2 = [1-x for x in data2]

popt2, pcov2 = curve_fit(function, xx22, data_can2, sigma=sd2, bounds=([0, 0, 1], [1., 5., 6.]))
perr2 = np.sqrt(np.diag(pcov2))  # standard deviation error - nie zawsze prawdziwe

print(popt2)
print(perr2)

residuals = data_can2 - function(xx22, *popt2)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_can2-np.mean(data_can2))**2)
r_squared = 1 - (ss_res / ss_tot)
print("r^2", r_squared)

text = f'{popt2[0]:.3f}' + " (" + f'{perr2[0]:.3f}' + ")\t" + f'{popt2[1]:.5f}' + " (" + f'{perr2[1]:.5f}' + ")\t" + f'{popt2[2]:.3f}' + " (" + f'{perr2[2]:.3f}' + ")\t" + f'{r_squared:.3f}' + "\n"
print(text)

ax[1].plot(xx2, [function(x, popt2[0], popt2[1], popt2[2]) for x in xx22], 'g--', label='fit: c=%.3f, a=%.3f, n=%.3f' % (popt2[0], popt2[1], popt2[2]))
ax[1].errorbar(xx2, data_can2, yerr=sd2, fmt='o', label="data")
ax[1].legend()
ax[1].set_xlabel("time (weeks)")
ax[1].set_ylabel("carcinoma probability genotype 1  - BRCA1 L63X/+")
plt.show()

# """  # avrami fit, paramteter c

# """


def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


def function_lin(x, a, n):
    return a + n * np.array(x)


x1 = [21, 38, 59, 76]
x1_scaled = [x*7/365 for x in x1]
data = [0.962, 0.92, 0.874, 0.794]
data_can = [1-x for x in data]
sd1 = [0.0377, 0.0545, 0.0685, 0.098]

x1_log = [np.log(x) for x in x1_scaled]
y1_log = [np.log(-np.log(1 - y)) for y in data_can]
sd1_log = [np.abs(sd1[i]/(np.log(1-data_can[i])*(1 - data_can[i]))) for i in range(len(sd1))]

pear = pearsonr([np.log(x) for x in x1_scaled], [np.log(-np.log(1 - y)) for y in data_can])
print(pear)

popt, pcov = curve_fit(function_lin, x1_log, y1_log, sigma=sd1_log)
perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe
print(popt, np.exp(popt[0]))

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(x1, [function(x, np.exp(popt[0]), popt[1]) for x in x1_scaled], 'g--', label='fit: a=%.5f, n=%.3f' % (np.exp(popt[0]), popt[1]))
# ax[0].scatter(x1_log, y1_log)
ax[0].errorbar(x1, data_can, yerr=sd1, fmt='o', label="data")
ax[0].legend()
ax[0].set_xlabel("time (weeks)")
ax[0].set_ylabel("carcinoma probability genotype 0 - BRCA1 +/+")


xx2 = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
xx22 = [x*7/365 for x in xx2]
data2 = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
sd2 = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
data_can2 = [1-x for x in data2]

x2_log = [np.log(x) for x in xx22]
y2_log = [np.log(-np.log(1 - y)) for y in data_can2]
sd2_log = [np.abs(sd2[i]/(np.log(1-data_can2[i])*(1 - data_can2[i]))) for i in range(len(sd2))]

pear = pearsonr([np.log(x) for x in xx22], [np.log(-np.log(1 - y)) for y in data_can2])
print(pear)

popt2, pcov2 = curve_fit(function_lin, x2_log, y2_log, sigma=sd2_log)
perr2 = np.sqrt(np.diag(pcov2))  # standard deviation error - nie zawsze prawdziwe

print(popt2, np.exp(popt2[0]))

ax[1].plot(xx2, [function(x, np.exp(popt2[0]), popt2[1]) for x in xx22], 'g--', label='fit: a=%.5f, n=%.3f' % (np.exp(popt2[0]), popt2[1]))
ax[1].errorbar(xx2, data_can2, yerr=sd2, fmt='o', label="data")
ax[1].legend()
ax[1].set_xlabel("time (weeks)")
ax[1].set_ylabel("carcinoma probability genotype 1  - BRCA1 L63X/+")
plt.show()

# """  # avrami fit c = 1

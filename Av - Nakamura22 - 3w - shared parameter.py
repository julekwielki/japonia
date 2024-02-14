import math

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


"""

def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


x1 = [21, 38, 59, 76]
x1_scaled = [x*7/365 for x in x1]
sd1 = [0.0377, 0.0545, 0.0685, 0.098]
data = [0.962, 0.92, 0.874, 0.794]
data_can = [1-x for x in data]

popt, pcov = curve_fit(function, x1_scaled, data_can, sigma=sd1, bounds=([0, 0], [5., 6.]))
perr = np.sqrt(np.diag(pcov))  # standard deviation error - nie zawsze prawdziwe


residuals = data_can - function(x1_scaled, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_can-np.mean(data_can))**2)
r_squared = 1 - (ss_res / ss_tot)

# print("genotype 0 - BRCA1 +/+")
# print(popt)
# print(perr)
# print("r^2", r_squared)

# text = "genotype 0 - BRCA1 +/+\na(u(a))\tn(u(n))\tr^2\n"
# text = text + f'{popt[0]:.3f}' + " (" + f'{perr[0]:.3f}' + ")\t" + f'{popt[1]:.5f}' + " (" + f'{perr[1]:.5f}' + ")\t" + f'{r_squared:.3f}' + "\n"

text = "genotype 0 - BRCA1 +/+\na\tu(a)\tn\tu(n)\tr^2\n"
text = text + f'{popt[0]:.3f}' + "\t" + f'{perr[0]:.3f}' + "\t" + f'{popt[1]:.5f}' + "\t" + f'{perr[1]:.5f}' + "\t" + f'{r_squared:.3f}' + "\n"
print(text)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(x1, [function(x, popt[0], popt[1]) for x in x1_scaled], 'g--', label='fit: a=%.5f, n=%.3f' % (popt[0], popt[1]))
ax[0].errorbar(x1, data_can, yerr=sd1, fmt='o', label="data")
ax[0].legend()
ax[0].set_xlabel("time (weeks)")
ax[0].set_ylabel("carcinoma probability genotype 0 - BRCA1 +/+")


xx2 = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
xx22 = [x*7/365 for x in xx2]
data2 = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
sd2 = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
data_can2 = [1-x for x in data2]

popt2, pcov2 = curve_fit(function, xx22, data_can2, sigma=sd2, bounds=([0, 0], [5., 8.]))
perr2 = np.sqrt(np.diag(pcov2))  # standard deviation error - nie zawsze prawdziwe


residuals = data_can2 - function(xx22, *popt2)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_can2-np.mean(data_can2))**2)
r_squared2 = 1 - (ss_res / ss_tot)

# print("genotype 1  - BRCA1 L63X/+")
# print(popt2)
# print(perr2)
# print("r^2", r_squared2)

# text = "genotype 1  - BRCA1 L63X/+\na(u(a))\tn(u(n))\tr^2\n"
# text = text + f'{popt2[0]:.3f}' + " (" + f'{perr2[0]:.3f}' + ")\t" + f'{popt2[1]:.5f}' + " (" + f'{perr2[1]:.5f}' + ")\t" + f'{r_squared2:.3f}' + "\n"

text = "genotype 1  - BRCA1 L63X/+\na\tu(a)\tn\tu(n)\tr^2\n"
text = text + f'{popt2[0]:.3f}' + "\t" + f'{perr2[0]:.3f}' + "\t" + f'{popt2[1]:.5f}' + "\t" + f'{perr2[1]:.5f}' + "\t" + f'{r_squared2:.3f}' + "\n"
print(text)

ax[1].plot(xx2, [function(x, popt2[0], popt2[1]) for x in xx22], 'g--', label='fit: a=%.3f, n=%.3f' % (popt2[0], popt2[1]))
ax[1].errorbar(xx2, data_can2, yerr=sd2, fmt='o', label="data")
ax[1].legend()
ax[1].set_xlabel("time (weeks)")
ax[1].set_ylabel("carcinoma probability genotype 1  - BRCA1 L63X/+")
plt.show()

# """

# dane do sklejenia

xx1 = [21, 38, 59, 76]
x1_scaled = [x*7/365 for x in xx1]
sd1 = [0.0377, 0.0545, 0.0685, 0.098]
data = [0.962, 0.92, 0.874, 0.794]
data_can = [1-x for x in data]


xx2 = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
x2_scaled = [x*7/365 for x in xx2]
data2 = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
sd2 = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
data_can2 = [1-x for x in data2]

# sklejanie

comboX = np.append(x1_scaled, x2_scaled)
comboY = np.append(data_can, data_can2)
comboSD = np.append(sd1, sd2)


def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


def fun1(x, a, n1, n2):
    return 1 - np.exp(-a * np.power(x, n1))


def fun2(x, a, n1, n2):
    return 1 - np.exp(-a * np.power(x, n2))


def comboFunc(comboData, a, n1, n2):
    extract1 = comboData[:len(data_can)]  # first data
    extract2 = comboData[len(data_can):]  # second data
    result1 = fun1(extract1, a, n1, n2)
    result2 = fun2(extract2, a, n1, n2)

    return np.append(result1, result2)


popt3, pcov3 = curve_fit(comboFunc, comboX, comboY, sigma=comboSD, bounds=([0, 0, 0], [5., 6., 8.]))
perr3 = np.sqrt(np.diag(pcov3))  # standard deviation error - nie zawsze prawdziwe

print(popt3)
print(perr3)


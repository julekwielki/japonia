import math

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

# """


def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


fig, ax = plt.subplots(1, 2, figsize=(12, 5))

x_mut = [42, 47, 57, 59, 63, 66, 80, 81, 88, 90, 91]
x_mut_scaled = [x*7/365 for x in x_mut]
data_mut = [0.957, 0.911, 0.854, 0.797, 0.74, 0.673, 0.561, 0.449, 0.336, 0.224, 0.112]
sd_mut = [0.0425, 0.0601, 0.0788, 0.0919, 0.1014, 0.1123, 0.1387, 0.1496, 0.1484, 0.1348, 0.1041]
data_mut_corect = [1-x for x in data_mut]

popt_mut, pcov_mut = curve_fit(function, x_mut_scaled, data_mut_corect, sigma=sd_mut, bounds=([0, 0], [5., 8.]))
perr_mut = np.sqrt(np.diag(pcov_mut))  # standard deviation error - nie zawsze prawdziwe


residuals = data_mut_corect - function(x_mut_scaled, *popt_mut)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_mut_corect-np.mean(data_mut_corect))**2)
r_squared_mut = 1 - (ss_res / ss_tot)

text = "genotype 1 - BRCA1 L63X/+\na\tu(a)\tn\tu(n)\tr^2\n"
text = text + f'{popt_mut[0]:.3f}' + "\t" + f'{perr_mut[0]:.3f}' + "\t" + f'{popt_mut[1]:.5f}' + "\t" + f'{perr_mut[1]:.5f}' + "\t" + f'{r_squared_mut:.3f}' + "\n"
print(text)

ax[1].errorbar(x_mut, data_mut_corect, yerr=sd_mut, fmt='o', label="data_norm mut")
ax[1].plot(x_mut, [function(x, popt_mut[0], popt_mut[1]) for x in x_mut_scaled], 'g--', label='fit: a=%.3f, n=%.3f' % (popt_mut[0], popt_mut[1]))
ax[1].legend()
ax[1].set_xlabel("time (weeks)")
ax[1].set_ylabel("carcinoma probability genotype 1 - BRCA1 L63X/+")


x_norm = [21, 38, 59, 76]
x_norm_scaled = [x*7/365 for x in x_norm]
sd_norm = [0.0377, 0.0545, 0.0685, 0.098]
data_norm = [0.962, 0.92, 0.874, 0.794]
data_norm_corect = [1-x for x in data_norm]

popt_norm, pcov_norm = curve_fit(function, x_norm_scaled, data_norm_corect, sigma=sd_norm, bounds=([0, 0], [5., 6.]))
perr_norm = np.sqrt(np.diag(pcov_norm))  # standard deviation error - nie zawsze prawdziwe


residuals = data_norm_corect - function(x_norm_scaled, *popt_norm)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_norm_corect-np.mean(data_norm_corect))**2)
r_squared_norm = 1 - (ss_res / ss_tot)

text = "genotype 0 - BRCA1 +/+\na\tu(a)\tn\tu(n)\tr^2\n"
text = text + f'{popt_norm[0]:.3f}' + "\t" + f'{perr_norm[0]:.3f}' + "\t" + f'{popt_norm[1]:.5f}' + "\t" + f'{perr_norm[1]:.5f}' + "\t" + f'{r_squared_norm:.3f}' + "\n"
print(text)

ax[0].errorbar(x_norm, data_norm_corect, yerr=sd_norm, fmt='o', label="data_norm norm")
ax[0].plot(x_norm, [function(x, popt_norm[0], popt_norm[1]) for x in x_norm_scaled], 'g--', label='fit: a=%.5f, n=%.3f' % (popt_norm[0], popt_norm[1]))
ax[0].legend()
ax[0].set_xlabel("time (weeks)")
ax[0].set_ylabel("carcinoma probability genotype 0 - BRCA1 +/+")
plt.show()

plt.errorbar(x_norm, data_norm_corect, yerr=sd_norm, fmt='o', label="data_norm norm")
plt.errorbar(x_mut, data_mut_corect, yerr=sd_mut, fmt='o', label="data_norm mut")
plt.show()


def function2(x, n):
    return 1 - np.exp(-popt_mut[0] * np.power(x, n))

popt_norm_z_mut, pcov_norm_z_mut = curve_fit(function2, x_norm_scaled, data_norm_corect, sigma=sd_norm, bounds=([0], [6.]))
perr_norm_z_mut = np.sqrt(np.diag(pcov_norm_z_mut))  # standard deviation error - nie zawsze prawdziwe


residuals = data_norm_corect - function2(x_norm_scaled, *popt_norm_z_mut)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((data_norm_corect-np.mean(data_norm_corect))**2)
r_squared_norm_z_mut = 1 - (ss_res / ss_tot)

text = "genotype 0 - BRCA1 +/+\nn\tu(n)\tr^2\n"
text = text + f'{popt_norm_z_mut[0]:.3f}' + "\t" + f'{perr_norm_z_mut[0]:.3f}' + "\t" + f'{r_squared_norm_z_mut:.3f}' + "\n"
print(text)

plt.errorbar(x_norm, data_norm_corect, yerr=sd_norm, fmt='o', label="data_norm norm")
plt.plot(x_norm, [function2(x, popt_norm_z_mut[0]) for x in x_norm_scaled], 'g--', label='fit: n=%.3f' % (popt_norm[1]))
plt.errorbar(x_mut, data_mut_corect, yerr=sd_mut, fmt='o', label="data_norm mut")
plt.plot(x_mut, [function(x, popt_mut[0], popt_mut[1]) for x in x_mut_scaled], 'g--', label='fit: a=%.3f, n=%.3f' % (popt_mut[0], popt_mut[1]))

plt.legend()
plt.xlabel("time (weeks)")
plt.ylabel("carcinoma probability")
plt.show()

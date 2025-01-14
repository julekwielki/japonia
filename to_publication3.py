
import Avrami_data
import matplotlib.pyplot as plt
import numpy as np
from Av_lin_fits_class import Fitted_one_log
import matplotlib.ticker


def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))

"""
data_all = []
text = ""
fig, ax = plt.subplots(2, 1, figsize=(8, 15))

time_all0 = Avrami_data.T_time_all0
time_all0.append(Avrami_data.T_20w_time70)
survival_all0 = Avrami_data.T_survival_all0
survival_all0.append(Avrami_data.T_20w_survival70)
std_all0 = Avrami_data.T_std_all0
std_all0.append(Avrami_data.T_20w_std70)

time_all1 = Avrami_data.T_time_all1
time_all1.append(Avrami_data.T_20w_time71)
survival_all1 = Avrami_data.T_survival_all1
survival_all1.append(Avrami_data.T_20w_survival71)
std_all1 = Avrami_data.T_std_all1
std_all1.append(Avrami_data.T_20w_std71)

name0 = ["non irradiated, virgin", "irradiated at 3 weeks, virgin", "irradiated at 7 weeks, virgin", "irradiated at 7 weeks, virgin, late"]
name1 = ["non irradiated, parous", "irradiated at 3 weeks, parous", "irradiated at 7 weeks, parous", "irradiated at 7 weeks, parous, late"]


col = ['blue', 'red', 'orange', 'plum', 'green', 'brown', 'olive', 'cyan']
col0 = ['blue', 'orange', 'green', 'olive']
col1 = ['red', 'plum', 'brown', 'cyan']

for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2'], name0[i], col0[i]]
    data_all.append(aaa)

    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i])

    ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2'], name1[i], col1[i]]
    data_all.append(aaa)

    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col1[i])

    ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i])

    text = text + name1[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

box = ax[0].get_position()
ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax[0].legend(loc='center left', bbox_to_anchor=(1, -0.1))
box = ax[1].get_position()
ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])

ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")
ax[1].set_xticks([0.5, 1, 1.6])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')
plt.show()

N_non_mut_x_years = Avrami_data.N_non_mut_x_years
N_non_mut_sd = Avrami_data.N_non_mut_sd
N_non_mut_canc = Avrami_data.N_non_mut_canc

A = Fitted_one_log(N_non_mut_x_years, N_non_mut_canc, N_non_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 normal", "teal"]
text = text + "BRCA1 normal" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

N_mut_x_years = Avrami_data.N_mut_x_years
N_mut_sd = Avrami_data.N_mut_sd
N_mut_canc = Avrami_data.N_mut_canc

A = Fitted_one_log(N_mut_x_years, N_mut_canc, N_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 mutation", "indigo"]
text = text + "BRCA1 mutation" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

x = []
y = []
for i in range(len(data_all)):
    plt.errorbar(data_all[i][0], data_all[i][2], label=data_all[i][5], xerr=data_all[i][1], yerr=data_all[i][3], fmt='o', ecolor=data_all[i][6], color=data_all[i][6])
plt.xlabel("parameter a")
plt.ylabel("parameter k")
plt.legend()
plt.show()

print(text) # """  # 1

"""
data_all = []
text = ""
fig, ax = plt.subplots(2, 1, figsize=(8, 15))

time_all0 = [Avrami_data.T_time70, Avrami_data.T_20w_time70]
survival_all0 = [Avrami_data.T_survival70, Avrami_data.T_20w_survival70]
std_all0 = [Avrami_data.T_std70, Avrami_data.T_20w_std70]

time_all1 = [Avrami_data.T_time71, Avrami_data.T_20w_time71]
survival_all1 = [Avrami_data.T_survival71, Avrami_data.T_20w_survival71]
std_all1 = [Avrami_data.T_std71, Avrami_data.T_20w_std71]

name0 = ["irradiated at 7 weeks, virgin", "irradiated at 7 weeks, virgin, late"]
name1 = ["irradiated at 7 weeks, parous", "irradiated at 7 weeks, parous, late"]


col = ['green', 'brown', 'olive', 'cyan']
col0 = ['green', 'olive']
col1 = ['brown', 'cyan']

for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2'], name0[i], col0[i]]
    data_all.append(aaa)

    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i])

    ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2'], name1[i], col1[i]]
    data_all.append(aaa)

    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col1[i])

    ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i])

    text = text + name1[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

box = ax[0].get_position()
ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax[0].legend(loc='center left', bbox_to_anchor=(1, -0.1))
box = ax[1].get_position()
ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])

ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")
ax[1].set_xticks([0.5, 1, 1.6])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')
plt.show()

N_non_mut_x_years = Avrami_data.N_non_mut_x_years
N_non_mut_sd = Avrami_data.N_non_mut_sd
N_non_mut_canc = Avrami_data.N_non_mut_canc

A = Fitted_one_log(N_non_mut_x_years, N_non_mut_canc, N_non_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 normal", "teal"]
text = text + "BRCA1 normal" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

N_mut_x_years = Avrami_data.N_mut_x_years
N_mut_sd = Avrami_data.N_mut_sd
N_mut_canc = Avrami_data.N_mut_canc

A = Fitted_one_log(N_mut_x_years, N_mut_canc, N_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 mutation", "indigo"]
text = text + "BRCA1 mutation" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

x = []
y = []
for i in range(len(data_all)):
    plt.errorbar(data_all[i][0], data_all[i][2], label=data_all[i][5], xerr=data_all[i][1], yerr=data_all[i][3], fmt='o', ecolor=data_all[i][6], color=data_all[i][6])
plt.xlabel("parameter a")
plt.ylabel("parameter k")
plt.legend()
plt.show()

print(text) # """  # 2

"""
data_all = []
text = ""
fig, ax = plt.subplots(2, 2, figsize=(8, 15))

time_all0 = [Avrami_data.T_time70]
survival_all0 = [Avrami_data.T_survival70]
std_all0 = [Avrami_data.T_std70]

time_all1 = [Avrami_data.T_time71]
survival_all1 = [Avrami_data.T_survival71]
std_all1 = [Avrami_data.T_std71]

name0 = ["irradiated at 7 weeks, virgin"]
name1 = ["irradiated at 7 weeks, parous"]


col = ['green', 'brown']
col0 = ['green']
col1 = ['brown']

for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2'], name0[i], col0[i]]
    data_all.append(aaa)

    ax[0][0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[0][0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i])

    ax[1][0].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[1][0].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2'], name1[i], col1[i]]
    data_all.append(aaa)

    ax[0][0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0][0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col1[i])

    ax[1][0].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1][0].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i])

    text = text + name1[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

time_all0 = [Avrami_data.T_20w_time70]
survival_all0 = [Avrami_data.T_20w_survival70]
std_all0 = [Avrami_data.T_20w_std70]

time_all1 = [Avrami_data.T_20w_time71]
survival_all1 = [Avrami_data.T_20w_survival71]
std_all1 = [Avrami_data.T_20w_std71]

name0 = ["irradiated at 7 weeks, virgin, late"]
name1 = ["irradiated at 7 weeks, parous, late"]


col = [ 'olive', 'cyan']
col0 = ['olive']
col1 = ['cyan']

for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2'], name0[i], col0[i]]
    data_all.append(aaa)

    ax[0][1].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[0][1].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i])

    ax[1][1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[1][1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2'], name1[i], col1[i]]
    data_all.append(aaa)

    ax[0][1].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0][1].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col1[i])

    ax[1][1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1][1].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i])

    text = text + name1[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"


box = ax[0][1].get_position()
ax[0][1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax[0][1].legend(loc='center left', bbox_to_anchor=(1, -0.1))
box = ax[1][1].get_position()
ax[1][1].set_position([box.x0, box.y0, box.width * 0.8, box.height])

ax[1][0].set_xscale('log')
ax[0][0].set_xlabel("time (years)")
ax[0][0].set_ylabel("Carcinoma incidence ")
ax[1][0].set_xlabel("time (years)")
ax[1][0].set_ylabel("linearised data")
ax[1][0].set_xticks([0.5, 1, 1.6])
ax[1][0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1][0].get_xaxis().set_tick_params(which='major')

ax[1][1].set_xscale('log')
ax[0][1].set_xlabel("time (years)")
ax[0][1].set_ylabel("Carcinoma incidence ")
ax[1][1].set_xlabel("time (years)")
ax[1][1].set_ylabel("linearised data")
ax[1][1].set_xticks([0.5, 1, 1.6])
ax[1][1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1][1].get_xaxis().set_tick_params(which='major')

plt.show()

N_non_mut_x_years = Avrami_data.N_non_mut_x_years
N_non_mut_sd = Avrami_data.N_non_mut_sd
N_non_mut_canc = Avrami_data.N_non_mut_canc

A = Fitted_one_log(N_non_mut_x_years, N_non_mut_canc, N_non_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 normal", "teal"]
text = text + "BRCA1 normal" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

N_mut_x_years = Avrami_data.N_mut_x_years
N_mut_sd = Avrami_data.N_mut_sd
N_mut_canc = Avrami_data.N_mut_canc

A = Fitted_one_log(N_mut_x_years, N_mut_canc, N_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 mutation", "indigo"]
text = text + "BRCA1 mutation" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

x = []
y = []
for i in range(len(data_all)):
    plt.errorbar(data_all[i][0], data_all[i][2], label=data_all[i][5], xerr=data_all[i][1], yerr=data_all[i][3], fmt='o', ecolor=data_all[i][6], color=data_all[i][6])
plt.xlabel("parameter a")
plt.ylabel("parameter k")
plt.legend()
plt.show()

print(text) # """  # 3

# """
data_all = []
text = ""
fig, ax = plt.subplots(2, 2, figsize=(8, 15))

time_all0 = [Avrami_data.T_time70]
survival_all0 = [Avrami_data.T_survival70]
std_all0 = [Avrami_data.T_std70]

time_all1 = [Avrami_data.T_20w_time70]
survival_all1 = [Avrami_data.T_20w_survival70]
std_all1 = [Avrami_data.T_20w_std70]

name0 = ["irradiated at 7 weeks, virgin"]
name1 = ["irradiated at 7 weeks, virgin, late"]


col = ['green', 'olive']
col0 = ['green']
col1 = ['olive']

for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2'], name0[i], col0[i]]
    data_all.append(aaa)

    ax[0][0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[0][0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i])

    ax[1][0].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[1][0].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2'], name1[i], col1[i]]
    data_all.append(aaa)

    ax[0][0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0][0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col1[i])

    ax[1][0].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1][0].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i])

    text = text + name1[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

time_all0 = [Avrami_data.T_time71]
survival_all0 = [Avrami_data.T_survival71]
std_all0 = [Avrami_data.T_std71]

time_all1 = [Avrami_data.T_20w_time71]
survival_all1 = [Avrami_data.T_20w_survival71]
std_all1 = [Avrami_data.T_20w_std71]

name0 = ["irradiated at 7 weeks, parous"]
name1 = ["irradiated at 7 weeks, parous, late"]


col = ['brown', 'cyan']
col0 = ['brown']
col1 = ['cyan']

for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2'], name0[i], col0[i]]
    data_all.append(aaa)

    ax[0][1].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[0][1].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i])

    ax[1][1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i], label=name0[i])
    ax[1][1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2'], name1[i], col1[i]]
    data_all.append(aaa)

    ax[0][1].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0][1].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col1[i])

    ax[1][1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1][1].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i])

    text = text + name1[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.4f}' + " " + f'{aaa[1]:.4f}' + "\t" + f'{aaa[2]:.2f}' + " " + f'{aaa[3]:.2f}' + "\t" + f'{aaa[4]:.3f}' + "\n"


box = ax[0][1].get_position()
ax[0][1].set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax[0][1].legend(loc='center left', bbox_to_anchor=(1, -0.1))
box = ax[1][1].get_position()
ax[1][1].set_position([box.x0, box.y0, box.width * 0.8, box.height])

ax[1][0].set_xscale('log')
ax[0][0].set_xlabel("time (years)")
ax[0][0].set_ylabel("Carcinoma incidence ")
ax[1][0].set_xlabel("time (years)")
ax[1][0].set_ylabel("linearised data")
ax[1][0].set_xticks([0.5, 1, 1.6])
ax[1][0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1][0].get_xaxis().set_tick_params(which='major')

ax[1][1].set_xscale('log')
ax[0][1].set_xlabel("time (years)")
ax[0][1].set_ylabel("Carcinoma incidence ")
ax[1][1].set_xlabel("time (years)")
ax[1][1].set_ylabel("linearised data")
ax[1][1].set_xticks([0.5, 1, 1.6])
ax[1][1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1][1].get_xaxis().set_tick_params(which='major')

plt.show()

N_non_mut_x_years = Avrami_data.N_non_mut_x_years
N_non_mut_sd = Avrami_data.N_non_mut_sd
N_non_mut_canc = Avrami_data.N_non_mut_canc

A = Fitted_one_log(N_non_mut_x_years, N_non_mut_canc, N_non_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 normal", "teal"]
text = text + "BRCA1 normal" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

N_mut_x_years = Avrami_data.N_mut_x_years
N_mut_sd = Avrami_data.N_mut_sd
N_mut_canc = Avrami_data.N_mut_canc

A = Fitted_one_log(N_mut_x_years, N_mut_canc, N_mut_sd)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "BRCA1 mutation", "indigo"]
text = text + "BRCA1 mutation" + "\na(u(a))\tn(u(n))\tr^2\n"
text = text + f'{aaa[0]:.4f}' + " (" + f'{aaa[1]:.4f}' + ")\t" + f'{aaa[2]:.2f}' + " (" + f'{aaa[3]:.2f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"
data_all.append(aaa)

x = []
y = []
for i in range(len(data_all)):
    plt.errorbar(data_all[i][0], data_all[i][2], label=data_all[i][5], xerr=data_all[i][1], yerr=data_all[i][3], fmt='o', ecolor=data_all[i][6], color=data_all[i][6])
plt.xlabel("parameter a")
plt.ylabel("parameter k")
plt.legend()
plt.show()

print(text) # """  # 4
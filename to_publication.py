
import Avrami_data
import matplotlib.pyplot as plt
import numpy as np
from Av_lin_fits_class import Fitted_one_log
import matplotlib.ticker

time_all0 = Avrami_data.T_20w_time_all0
survival_all0 = Avrami_data.T_20w_survival_all0
std_all0 = Avrami_data.T_20w_std_all0

time_all1 = Avrami_data.T_20w_time_all1
survival_all1 = Avrami_data.T_20w_survival_all1
std_all1 = Avrami_data.T_20w_std_all1


# """
time_all0 = Avrami_data.T_time_all0
survival_all0 = Avrami_data.T_survival_all0
std_all0 = Avrami_data.T_std_all0

time_all1 = Avrami_data.T_time_all1
survival_all1 = Avrami_data.T_survival_all1
std_all1 = Avrami_data.T_std_all1
# """
name0 = ["t0 p0", "t3 p0", "t7 p0"]
name0 = ["non irradiated, virgin", "irradiated at 3 weeks, virgin", "irradiated at 7 weeks, virgin"]
name1 = ["t0 p1", "t3 p1", "t7 p1"]
name1 = ["non irradiated, parous", "irradiated at 3 weeks, parous", "irradiated at 7 weeks, parous"]
namename = ["non irradiated, virgin", "non irradiated, parous", "irradiated at 3 weeks, virgin", "irradiated at 3 weeks, parous", "irradiated at 7 weeks, virgin", "irradiated at 7 weeks, parous"]
col = ['blue', 'red', 'orange', 'plum', 'green', 'brown']
col0 = ['blue', 'orange', 'green']
col1 = ['red', 'plum', 'brown']

data_all = []
text = ""

fig, ax = plt.subplots(2, 1, figsize=(8, 15))


def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


for i in range(len(time_all0)):
    xx2 = time_all0[i]
    xx = [x*7/365 for x in xx2]
    sd = std_all0[i]
    data = [1-x for x in survival_all0[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0]*np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1], dane0['r2']]

    data_all.append(aaa)

    # """
    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col0[i], color=col0[i])
    ax[0].scatter(xx, data, label=name0[i], color=col0[i])
    ax[0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx], '--', color=col0[i]) # , label='fit: n=%.1f' % (popt[1]), color=col[i])"""

    ax[1].scatter(xx, A.y, label=name0[i], color=col0[i])
    ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col0[i], color=col0[i])
    ax[1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x],
               '--', color=col0[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    """  xx = time_all1[i]
    # xx2 = [x*7/365 for x in xx]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    data_all.append(dane1)

    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i])
    ax[0].scatter(xx, data, label=name1[i], color=col1[i])
    ax[0].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x], '--', color=col1[i]) # , label='fit: n=%.1f' % (popt[1]), color=col[i])

    ax[1].scatter(xx, [np.log(-np.log(1 - y)) for y in data], label=name1[i], color=col1[i])
    ax[1].errorbar(xx, [np.log(-np.log(1 - y)) for y in data],
                   [np.abs(sd[i] / (np.log(1 - data[i]) * (1 - data[i]))) for i in range(len(sd))], fmt='o',
                   ecolor=col1[i], color=col1[i])"""

    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.6f}' + " (" + f'{aaa[1]:.6f}' + ")\t" + f'{aaa[2]:.5f}' + " (" + f'{aaa[3]:.5f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"

    xx2 = time_all1[i]
    xx = [x * 7 / 365 for x in xx2]
    sd = std_all1[i]
    data = [1 - x for x in survival_all1[i]]

    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane1 = A.toasty()
    aaa = [np.exp(dane1["param"][0]), dane1['err'][0] * np.exp(dane1["param"][0]), dane1["param"][1], dane1['err'][1],
           dane1['r2']]

    data_all.append(aaa)

    # """
    ax[0].errorbar(xx, data, yerr=sd, fmt='o', ecolor=col1[i], color=col1[i])
    ax[0].scatter(xx, data, label=name1[i], color=col1[i])
    ax[0].plot(xx, [function(x, aaa[0], aaa[2]) for x in xx],
               '--', color=col1[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])"""

    ax[1].scatter(xx, A.y, label=name1[i], color=col1[i])
    ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i])
    ax[1].plot(xx, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in A.x],
               '--', color=col1[i])  # , label='fit: n=%.1f' % (popt[1]), color=col[i])


    text = text + name0[i] + "\na(u(a))\tn(u(n))\tr^2\n"
    text = text + f'{aaa[0]:.6f}' + " (" + f'{aaa[1]:.6f}' + ")\t" + f'{aaa[2]:.5f}' + " (" + f'{aaa[3]:.5f}' + ")\t" + f'{aaa[4]:.3f}' + "\n"


print(text)

# ax[0].legend(loc=2)
# ax[1].legend(loc=2)

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
ax[1].get_xaxis().set_tick_params(which='major')
plt.show()

x = []
y = []
for i in range(len(data_all)):
    print(i, data_all[i][0], data_all[i][2], namename[i], col[i])
    plt.errorbar(data_all[i][0], data_all[i][2], label=namename[i], xerr=data_all[i][1], yerr=data_all[i][3], fmt='o', ecolor=col[i], color=col[i])
plt.legend()
plt.show()


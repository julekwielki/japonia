from scipy.stats import ks_2samp

import Avrami_data
from Av_lin_fits_class import Fitted
import matplotlib.pyplot as plt
import numpy as np


"""x1 = Avrami_data.N_non_mut_x_years
x2 = Avrami_data.N_mut_x_years

y1 = Avrami_data.N_non_mut_canc
y2 = Avrami_data.N_mut_canc

sd1 = Avrami_data.N_non_mut_sd
sd2 = Avrami_data.N_mut_sd

A = Fitted(x1, x2, y1, y2, sd1, sd2, [0, 0, 0], [1, 10, 10])
A.final_fit_sd()
A.plot_data_sd()
a, b, c = A.toasty()
print(a, b, c)"""

"""
time_all0 = Avrami_data.T_time_all0
survival_all0 = Avrami_data.T_survival_all0
std_all0 = Avrami_data.T_std_all0

time_all1 = Avrami_data.T_time_all1
survival_all1 = Avrami_data.T_survival_all1
std_all1 = Avrami_data.T_std_all1
# """

# """
time_all0 = Avrami_data.T_20w_time_all0
survival_all0 = Avrami_data.T_20w_survival_all0
std_all0 = Avrami_data.T_20w_std_all0

time_all1 = Avrami_data.T_20w_time_all1
survival_all1 = Avrami_data.T_20w_survival_all1
std_all1 = Avrami_data.T_20w_std_all1
# """

data1 = []
data2 = []
names = ["t0 p0", "t0 p1", "t3 p1", "t3 p0", "t7 p1", "t7 p0"]
names = ["t0 p0", "t3 p0", "t7 p0", "t0 p1", "t3 p1", "t7 p1"]
col = ['blue', 'orange', 'green', 'red', 'plum', 'brown']
# fig, ax = plt.subplots(1, 2, figsize=(12, 8))
fig, ax = plt.subplots(2, 1, figsize=(8, 15))
# """
for i in range(len(time_all0)):
    x1 = [x * 7/365 for x in time_all0[i]]
    x2 = [x * 7/365 for x in time_all1[i]]

    y1 = [1-x for x in survival_all0[i]]
    y2 = [1-x for x in survival_all1[i]]

    sd1 = std_all0[i]
    sd2 = std_all1[i]

    # """
    print(names[i])

    A = Fitted(x1, y1, sd1)
    A.fit_both_sd()
    # A.plot_data_sd(names[i])
    dane1, dane2 = A.toasty()
    data1.append(dane1)
    data2.append(dane2)

    # """
    ax[0].scatter(x1, y1, label=names[i], color=col[i])
    ax[0].errorbar(x1, y1, sd1, fmt='o', ecolor=col[i], color=col[i])  # """

    """
    ax[0].plot(x1, [A.function(x, dane1["param"][0], dane1["param"][1]) for x in x1], '--',
               label='fit norm: a=%.3f, n=%.1f' % (dane1["param"][0], dane1["param"][1]), color=col[i]) # """
    """
    ax[0].plot(x1, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x1], '--',
               label='fit ln: a=%.3f, n=%.1f' % (dane2["param"][0], dane2["param"][1]), color=col[i])  # """
    # """
    ax[0].plot(x1, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x1], '--', color=col[i]) # """

    # """
    ax[1].scatter(x1, [np.log(-np.log(1 - y)) for y in y1], label=names[i], color=col[i])
    ax[1].errorbar(x1, [np.log(-np.log(1 - y)) for y in y1],
                   [np.abs(sd1[i] / (np.log(1 - y1[i]) * (1 - y1[i]))) for i in range(len(sd1))], fmt='o',
                   ecolor=col[i], color=col[i])  # """
    """
    ax[1].plot(x1, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in x1], '--',
               label='fit norm: a=%.3f, n=%.1f' % (dane1["param"][0], dane1["param"][1]), color=col[i])  # """
    """
    ax[1].plot(x1, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x1], '--',
               label='fit ln: a=%.3f, n=%.1f' % (dane2["param"][0], dane2["param"][1]), color=col[i])  # """

    # """
    ax[1].plot(x1, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x1], '--', color=col[i])  # """

    """
    fig, ax = plt.subplots(2, 2, figsize=(12, 12))

    ax[0][0].scatter(x1, y1, label=names[i], color='orange')
    ax[0][0].errorbar(x1, y1, sd1, fmt='o', ecolor='orange', color='orange')
    ax[0][0].plot(x1, [A.function(x, dane1["param"][0], dane1["param"][1]) for x in x1], '--',
                  label='fit norm: a=%.5f, n=%.3f' % (dane1["param"][0], dane1["param"][1]), color='orange')
    ax[0][0].plot(x1, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x1], '--',
                  label='fit ln: a=%.5f, n=%.3f' % (dane2["param"][0], dane2["param"][1]), color='blue')
    ax[0][0].legend()
    ax[0][0].set_ylabel("carcinoma probability")
    ax[0][0].set_xlabel("time (years)")

    ax[0][1].scatter(x1, [np.log(-np.log(1 - y)) for y in y1], label=names[i], color='blue')
    ax[0][1].errorbar(x1, [np.log(-np.log(1 - y)) for y in y1], [np.abs(sd1[i]/(np.log(1-y1[i])*(1 - y1[i]))) for i in range(len(sd1))], fmt='o', ecolor='blue', color='blue')
    ax[0][1].plot(x1, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x1], '--',
                  label='fit ln: a=%.5f, n=%.3f' % (dane2["param"][0], dane2["param"][1]), color='blue')

    ax[0][1].plot(x1, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in x1], '--',
                  label='fit norm: a=%.5f, n=%.3f' % (dane1["param"][0], dane1["param"][1]), color='orange')

    ax[0][1].set_xscale('log')
    ax[0][1].legend()
    ax[0][1].set_xlabel("time (years)")
    ax[0][1].set_ylabel("zlogarytmowane")  # """

    """
    print(names[i])
    print("a\tu(a)\tn\tu(n)\tr2")
    print(dane1["param"][0], dane1["err"][0], dane1["param"][1], dane1["err"][1], dane1["r2"])
    print(dane2["param"][0], dane2["err"][0], dane2["param"][1], dane2["err"][1], dane2["r2"])  # """

    print(names[i+3])

    B = Fitted(x2, y2, sd2)
    B.fit_both_sd()
    # B.plot_data_sd(names[i+3])
    dane1, dane2 = B.toasty()
    data1.append(dane1)
    data2.append(dane2)

    # """
    ax[0].errorbar(x2, y2, sd2, fmt='o', ecolor=col[i+3], color=col[i+3])
    ax[0].scatter(x2, y2, label=names[i+3], color=col[i+3])  # """

    """
    ax[0].plot(x2, [A.function(x, dane1["param"][0], dane1["param"][1]) for x in x2], '--',
               label='fit norm: a=%.3f, n=%.1f' % (dane1["param"][0], dane1["param"][1]), color=col[i+3])  # """
    """
    ax[0].plot(x2, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--',
               label='fit ln: a=%.3f, n=%.1f' % (dane2["param"][0], dane2["param"][1]), color=col[i + 3])  # """
    # """
    ax[0].plot(x2, [A.function(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--', color=col[i + 3])  # """

    # """
    ax[1].scatter(x2, [np.log(-np.log(1 - y)) for y in y2], label=names[i+3], color=col[i+3])
    ax[1].errorbar(x2, [np.log(-np.log(1 - y)) for y in y2],
                   [np.abs(sd2[i] / (np.log(1 - y2[i]) * (1 - y2[i]))) for i in range(len(sd2))], fmt='o',
                   ecolor=col[i+3], color=col[i+3])  # """
    """
    ax[1].plot(x2, [A.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in x2], '--',
               label='fit norm: a=%.3f, n=%.1f' % (dane1["param"][0], dane1["param"][1]), color=col[i+3])  # """
    """
    ax[1].plot(x2, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--',
               label='fit ln: a=%.3f, n=%.1f' % (dane2["param"][0], dane2["param"][1]), color=col[i+3]) # """
    # """
    ax[1].plot(x2, [A.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--', color=col[i + 3])  # """

    """
    print(names[i+3])
    print(dane1["param"][0], dane1["err"][0], dane1["param"][1], dane1["err"][1], dane1["r2"])
    print(dane2["param"][0], dane2["err"][0], dane2["param"][1], dane2["err"][1], dane2["r2"]) # """

    """
    ax[1][0].scatter(x2, y2, label=names[i+3], color='orange')
    ax[1][0].errorbar(x2, y2, sd2, fmt='o', ecolor='orange', color='orange')
    ax[1][0].plot(x2, [B.function(x, dane1["param"][0], dane1["param"][1]) for x in x2], '--',
                  label='fit norm: a=%.5f, n=%.3f' % (dane1["param"][0], dane1["param"][1]), color='orange')
    ax[1][0].plot(x2, [B.function(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--',
                  label='fit ln: a=%.5f, n=%.3f' % (dane2["param"][0], dane2["param"][1]), color='blue')
    ax[1][0].legend()
    ax[1][0].set_ylabel("carcinoma probability")
    ax[1][0].set_xlabel("time (years)")

    ax[1][1].scatter(x2, [np.log(-np.log(1 - y)) for y in y2], label=names[i+3], color='blue')
    ax[1][1].errorbar(x2, [np.log(-np.log(1 - y)) for y in y2],
                      [np.abs(sd2[i] / (np.log(1 - y2[i]) * (1 - y2[i]))) for i in range(len(sd2))], fmt='o',
                      ecolor='blue', color='blue')
    ax[1][1].plot(x2, [B.function_ln(x, dane2["param"][0], dane2["param"][1]) for x in x2], '--',
                  label='fit ln: a=%.5f, n=%.3f' % (dane2["param"][0], dane2["param"][1]), color='blue')

    ax[1][1].plot(x2, [B.function_ln(x, dane1["param"][0], dane1["param"][1]) for x in x2], '--',
                  label='fit norm: a=%.5f, n=%.3f' % (dane1["param"][0], dane1["param"][1]), color='orange')

    ax[1][1].set_xscale('log')
    ax[1][1].legend()
    ax[1][1].set_xlabel("time (years)")
    ax[1][1].set_ylabel("zlogarytmowane")
    plt.savefig(names[i] + " sd.png")
    plt.show()  # """
    plt.show()

# """
ax[0].legend(loc=2)
ax[1].legend(loc=2)
ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")
# ax[0].plot(x2, [A.function(x, 0.8133673, 1.6234204) for x in x2], '--', color='black')
# ax[1].plot(x2, [A.function_ln(x, 0.8133673, 1.6234204) for x in x2], '--', color='black')
# fig.suptitle("dopasowanie do danych zlogarytmowanych, dane Tabatake po 20 tygodniu, narysowane dla krzywej sigmoidalnej i dla danych zlogarytmowanych")

# plt.show()# """

for i in data2:
    print('a = %.3f (%.3f), k = %.3f (%.3f)' % (i["param"][0], i["err"][0], i["param"][1], i["err"][1]))


"""
aa_norm = []
nn_norm = []
aa_log = []
nn_log = []
names = ["t0 p0", "t0 p1", "t3 p0", "t3 p1", "t7 p0", "t7 p1"]
col = ['blue', 'red', 'orange', 'plum', 'green', 'brown']

for i in range(len(data1)):
    a_norm = data1[i]["param"][0]
    a_log = data2[i]["param"][0]

    n_norm = data1[i]["param"][1]
    n_log = data2[i]["param"][1]
    aa_norm.append(a_norm)
    nn_norm.append(n_norm)
    aa_log.append(a_log)
    nn_log.append(n_log)

    print(a_norm, a_log, n_norm, n_log, names[i])
    
    if names[i] != "t0 p1":
        ax[2].errorbar(a_log, n_log, xerr=data2[i]["err"][0], yerr=data2[i]["err"][1], label=names[i], ecolor=col[i], color=col[i])  
    
    # if names[i] != "t0 p1":
    #     plt.errorbar(a_norm, n_norm, xerr=data1[i]["err"][0], yerr=data1[i]["err"][1], label=names[i], ecolor=col[i], color=col[i])  

# plt.scatter(aa_norm, nn_norm)
ax[0].legend()
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("prawdopodobieństwo nowotworu")
ax[1].legend()
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("zlogarytmowane")
ax[1].set_xscale('log')
ax[2].legend()
ax[2].set_xlabel("parametr a")
ax[2].set_ylabel("parametr n")
plt.show()
# """

"""

aa_norm = []
nn_norm = []
aa_log = []
nn_log = []
names = ["t0 p0", "t0 p1", "t3 p0", "t3 p1", "t7 p0", "t7 p1"]
col = ['blue', 'red', 'orange', 'plum', 'green', 'brown']

for i in range(len(data1)):
    a_norm = data1[i]["param"][0]
    a_log = data2[i]["param"][0]

    n_norm = data1[i]["param"][1]
    n_log = data2[i]["param"][1]
    aa_norm.append(a_norm)
    nn_norm.append(n_norm)
    aa_log.append(a_log)
    nn_log.append(n_log)

    print(a_norm, a_log, n_norm, n_log, names[i])

    if names[i] != "t0 p1":
        plt.errorbar(a_log, n_log, xerr=data2[i]["err"][0], yerr=data2[i]["err"][1], label=names[i], ecolor=col[i],
                       color=col[i])

        # if names[i] != "t0 p1":
    #     plt.errorbar(a_norm, n_norm, xerr=data1[i]["err"][0], yerr=data1[i]["err"][1], label=names[i], ecolor=col[i], color=col[i])

plt.legend()
plt.xlabel("parametr a")
plt.ylabel("parametr n")
plt.show()

# """
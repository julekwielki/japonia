import math

import matplotlib.ticker
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from Av_lin_fits_class import Fitted_one_log
from matplotlib import animation
from matplotlib.animation import PillowWriter
from matplotlib.widgets import Slider,  Button

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
col1 = ['blue', 'red', 'orange', 'plum', 'green', 'brown']
name1 = ['control', '3mGy', '6mGy', '12mGy', '24mGy', '60mGy']




def function(x, a, b, c):
    return np.exp(- a * x**2 * np.exp(-b * x)) + c * x


def function_av(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


def function_sig(x, b, c):
    return 1 / (1 + np.exp(-b * (x - c)))


def animate_one_set(datax, datay, minim=20, xlab="dose rate", ylab="relative risk"):  # slider datay (datax)

    fig = plt.figure(3)
    ax1 = fig.add_subplot()
    ax1.plot(datax, datay[minim])

    plt.subplots_adjust(bottom=0.25)  # making space for slider
    ax1.set_title("week: " + str(minim))
    lim = [0, max(max(p) for p in datay[minim:])*1.1]
    ax1.set_ylim(lim)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)

    axt = plt.axes([0.1, 0.07, 0.8, 0.02])  # współrzędne slidera
    st = Slider(axt, 't', minim, len(datay), valstep=1)

    def update2(val):  #
        ax1.clear()
        ax1.plot(datax, datay[st.val - 1])
        ax1.set_ylim(lim)
        ax1.set_title("week: " + str(st.val))
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)

    st.on_changed(update2)
    plt.show()


def animate_two_sets(datax, datay, datax2, datay2, minim=20, xlab="dose rate", ylab="relative risk"):
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.plot(datax, datay[minim])
    ax1.plot(datax2, datay2[minim], 'g--')

    ax1.set_title("week: " + str(minim))
    lim = [0, max(max(p) for p in datay[minim:])*1.1]
    ax1.set_ylim(lim)
    ax1.plot(datax, np.ones_like(datax))
    ax1.set_title(minim)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)

    plt.subplots_adjust(bottom=0.25)
    axt = plt.axes([0.1, 0.07, 0.8, 0.02])
    st = Slider(axt, 't', minim, len(datay), valstep=1)

    def update2(val):
        ax1.clear()
        ax1.set_ylim(lim)
        ax1.set_title("week: " + str(st.val))
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.plot(datax, datay[st.val - 1])
        ax1.plot(datax, np.ones_like(datax))
        ax1.plot(datax2, datay2[st.val - 1], 'g--')

    st.on_changed(update2)
    plt.show()



sig_fun_params = [[], []]
av_fun_params = [[], []]

# fig, ax = plt.subplots(2, 1, figsize=(8, 7))

for i in range(len(dose_rates)):
    data = [1 - x for x in YY[i]]
    xx = XX[i]
    sd = SD[i]
    A = Fitted_one_log(xx, data, sd)
    A.fit_both_sd()
    dane0 = A.toasty()
    aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
           dane0['r2']]
    av_fun_params[0].append(aaa[0])
    av_fun_params[1].append(aaa[2])

    popt1, pcov1 = curve_fit(function_sig, XX[i], YY[i], bounds=([-1, 0], [0., 150.]), sigma=SD[i])
    sig_fun_params[0].append(popt1[0])
    sig_fun_params[1].append(popt1[1])

    residuals = YY[i] - function_sig(XX[i], *popt1)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((YY[i] - np.mean(YY[i])) ** 2)
    r_squared1 = 1 - (ss_res / ss_tot)


    """ax[0].errorbar(xx, YY[i], yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[0].plot(range(100), [1 - function_av(x, aaa[0], aaa[2]) for x in range(100)], '--', color=col1[i])

    ax[1].errorbar(xx, YY[i], yerr=sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    ax[1].plot(range(100), [function_sig(a, popt1[0], popt1[1]) for a in range(100)], '--', color=col1[i])"""

    # ax[1].errorbar(xx, A.y, A.sd, fmt='o', ecolor=col1[i], color=col1[i], label=name1[i])
    # ax[1].plot(xx, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x], '--', color=col1[i])
# ax[1].set_xscale('log')

"""box = ax[0].get_position()
ax[0].set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax[0].legend(loc='center left', bbox_to_anchor=(1.05, -0.1))
box = ax[1].get_position()
ax[1].set_position([box.x0, box.y0, box.width * 0.8, box.height])

ax[1].set_xticks([20, 30, 40, 50, 60, 80])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')

plt.show()"""

survival_from_function_sig = []  # survival = function_sig(time, b, c) - for given dose rat
cancer_risk_from_function_sig = []  # risk = 1 - survival
for x in range(100):
    a = []
    b = []
    for y in range(len(sig_fun_params[0])):  # for every dose rate
        i = function_sig(x, sig_fun_params[0][y], sig_fun_params[1][y])
        a.append(1 - i)
        b.append(i)
    survival_from_function_sig.append(a)
    cancer_risk_from_function_sig.append(b)


survival_from_function_av = []  # survival = function_sig(time, b, c) - for given dose rat
cancer_risk_from_function_av = []  # risk = 1 - survival
for x in range(1, 101, 1):
    a = []
    b = []
    for y in range(len(av_fun_params[0])):  # for every dose rate
        i = function_av(x, av_fun_params[0][y], av_fun_params[1][y])
        a.append(1 - i)
        b.append(i)
    survival_from_function_av.append(b)
    cancer_risk_from_function_av.append(a)

only_doses = [3, 6, 12, 24, 60]
do_wizki_RR = []  # cancer risk (CR) for a dose / CR for 0 Gy

for x in range(len(survival_from_function_sig)):
    a = []
    for y in range(len(only_doses)):
        a.append(survival_from_function_sig[x][y+1]/survival_from_function_sig[x][0])
    do_wizki_RR.append(a)  # RR for given week for each dose rate

#animate_one_set(only_doses, do_wizki_RR, 50)

abc = [[], [], []]
mi, ma = 30, 100
for nr in range(mi, ma):
    popt1, pcov1 = curve_fit(function, only_doses, do_wizki_RR[nr], bounds=([-0, 0, 0], [2., 2., 1]))
    abc[0].append(popt1[0])
    abc[1].append(popt1[1])
    abc[2].append(popt1[2])
    # print(nr, popt1[0], popt1[1], popt1[2])
    # plt.plot(range(60), [function(x, popt1[0], popt1[1], popt1[2]) for x in range(60)], label='fit: a=%.3f, b=%.3f, c=%.3f' % tuple(popt1))
    # plt.scatter(only_doses, do_wizki_RR[nr], label="data")
    # plt.legend()
    # plt.show()
# """
plt.title("exp(- a * x^2 * exp(-b * x)) + c * x")
plt.scatter(range(mi, ma), abc[0], label="a")
plt.scatter(range(mi, ma), abc[1], label="b")
plt.scatter(range(mi, ma), abc[2], label="c")
plt.legend()
plt.show()# """

dat = []
for nr in range(100):
    popt1, pcov1 = curve_fit(function, only_doses, do_wizki_RR[nr], bounds=([-0, 0, 0], [2., 2., 1]))
    # print(nr, popt1[0], popt1[1], popt1[2])
    a = [function(x, popt1[0], popt1[1], popt1[2]) for x in range(60)]
    dat.append(a)

animate_two_sets(only_doses, do_wizki_RR, range(60), dat, 30)

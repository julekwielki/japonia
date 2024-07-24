import Avrami_data
import matplotlib.pyplot as plt
import numpy as np
from Av_lin_fits_class import Fitted_one_log
import matplotlib.ticker


def function(x, a, n):
    return 1 - np.exp(-a * np.power(x, n))


fig, ax = plt.subplots(2, 1, figsize=(8, 15))

T_time71 = Avrami_data.T_time71
T_survival71 = Avrami_data.T_survival71
T_std71 = Avrami_data.T_std71

T_time70 = Avrami_data.T_time70
T_survival70 = Avrami_data.T_survival70
T_std70 = Avrami_data.T_std70

time_w = T_time70 + T_time71
se = T_std70 + T_std71
surv = T_survival70 + T_survival71

time = [x*7/365 for x in time_w]
data = [1-x for x in surv]
time, data, se = zip(*sorted(zip(time, data, se)))

A = Fitted_one_log(time, data, se)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "all_data", "teal"]

# ax[0].scatter(N_non_mut_x_years, N_non_mut_canc, label="BRCA1 normal", color="teal")
ax[0].errorbar(time, data, yerr=se, fmt='o', ecolor="blue", color="blue", label="irradiated at 7 weeks, all")
ax[0].plot(time, [function(x, aaa[0], aaa[2]) for x in time], '--', color="blue")

# ax[1].scatter(N_non_mut_x_years, A.y, label="BRCA1 normal", color="teal")
ax[1].errorbar(time, A.y, A.sd, fmt='o', ecolor="blue", color="blue", label="irradiated at 7 weeks, all")
ax[1].plot(time, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x], '--', color="blue")

ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")

ax[1].set_xticks([0.5, 1, 1.5])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')
plt.show()

fig, ax = plt.subplots(2, 1, figsize=(8, 15))

T_time71 = Avrami_data.T_20w_time71
T_survival71 = Avrami_data.T_20w_survival71
T_std71 = Avrami_data.T_20w_std71

T_time70 = Avrami_data.T_20w_time70
T_survival70 = Avrami_data.T_20w_survival70
T_std70 = Avrami_data.T_20w_std70

time_w = T_time70 + T_time71
se = T_std70 + T_std71
surv = T_survival70 + T_survival71

time = [x*7/365 for x in time_w]
data = [1-x for x in surv]
time, data, se = zip(*sorted(zip(time, data, se)))

A = Fitted_one_log(time, data, se)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "all_data", "teal"]

# ax[0].scatter(N_non_mut_x_years, N_non_mut_canc, label="BRCA1 normal", color="teal")
ax[0].errorbar(time, data, yerr=se, fmt='o', ecolor="blue", color="blue", label="irradiated at 7 weeks, late")
ax[0].plot(time, [function(x, aaa[0], aaa[2]) for x in time], '--', color="blue")

# ax[1].scatter(N_non_mut_x_years, A.y, label="BRCA1 normal", color="teal")
ax[1].errorbar(time, A.y, A.sd, fmt='o', ecolor="blue", color="blue", label="irradiated at 7 weeks, late")
ax[1].plot(time, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x], '--', color="blue")

ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")

ax[1].set_xticks([0.5, 1, 1.5])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')
plt.show()

fig, ax = plt.subplots(2, 1, figsize=(8, 15))

T_time71 = Avrami_data.T_time71
T_survival71 = Avrami_data.T_survival71
T_std71 = Avrami_data.T_std71

T_time70 = Avrami_data.T_time70
T_survival70 = Avrami_data.T_survival70
T_std70 = Avrami_data.T_std70

time_w = T_time70 + T_time71
se = T_std70 + T_std71
surv = T_survival70 + T_survival71

time = [x*7/365 for x in time_w]
data = [1-x for x in surv]
time, data, se = zip(*sorted(zip(time, data, se)))

A = Fitted_one_log(time, data, se)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "all_data", "teal"]

# ax[0].scatter(N_non_mut_x_years, N_non_mut_canc, label="BRCA1 normal", color="teal")
ax[0].errorbar(time, data, yerr=se, fmt='o', ecolor="green", color="green", label="irradiated at 7 weeks, all")
ax[0].plot(time, [function(x, aaa[0], aaa[2]) for x in time], '--', color="green")

# ax[1].scatter(N_non_mut_x_years, A.y, label="BRCA1 normal", color="teal")
ax[1].errorbar(time, A.y, A.sd, fmt='o', ecolor="green", color="green", label="irradiated at 7 weeks, all")
ax[1].plot(time, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x], '--', color="green")

ax[1].set_xscale('log')
ax[0].set_xlabel("time (years)")
ax[0].set_ylabel("Carcinoma incidence ")
ax[1].set_xlabel("time (years)")
ax[1].set_ylabel("linearised data")

ax[1].set_xticks([0.5, 1, 1.5])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')

# """


T_time71 = Avrami_data.T_20w_time71
T_survival71 = Avrami_data.T_20w_survival71
T_std71 = Avrami_data.T_20w_std71

T_time70 = Avrami_data.T_20w_time70
T_survival70 = Avrami_data.T_20w_survival70
T_std70 = Avrami_data.T_20w_std70

time_w = T_time70 + T_time71
se = T_std70 + T_std71
surv = T_survival70 + T_survival71

time = [x*7/365 for x in time_w]
data = [1-x for x in surv]
time, data, se = zip(*sorted(zip(time, data, se)))

A = Fitted_one_log(time, data, se)
A.fit_both_sd()
dane0 = A.toasty()
aaa = [np.exp(dane0["param"][0]), dane0['err'][0] * np.exp(dane0["param"][0]), dane0["param"][1], dane0['err'][1],
       dane0['r2'], "all_data", "teal"]

# ax[0].scatter(N_non_mut_x_years, N_non_mut_canc, label="BRCA1 normal", color="teal")
ax[0].errorbar(time, data, yerr=se, fmt='o', ecolor="brown", color="brown", label="irradiated at 7 weeks, late")
ax[0].plot(time, [function(x, aaa[0], aaa[2]) for x in time], '--', color="brown")

# ax[1].scatter(N_non_mut_x_years, A.y, label="BRCA1 normal", color="teal")
ax[1].errorbar(time, A.y, A.sd, fmt='o', ecolor="brown", color="brown", label="irradiated at 7 weeks, late")
ax[1].plot(time, [A.function_ln(x, dane0["param"][0], dane0["param"][1]) for x in A.x], '--', color="brown")

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

ax[1].set_xticks([0.5, 1, 1.5])
ax[1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax[1].get_xaxis().set_tick_params(which='major')
plt.show()

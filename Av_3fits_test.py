import Avrami_data
from Avrami_fits_class import Fitted
import matplotlib.pyplot as plt

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


x1 = [x * 7/365 for x in Avrami_data.T_20w_time70]
x2 = [x * 7/365 for x in Avrami_data.T_20w_time71]

y1 = [1-x for x in Avrami_data.T_20w_survival70]
y2 = [1-x for x in Avrami_data.T_20w_survival71]

sd1 = Avrami_data.T_20w_std70
sd2 = Avrami_data.T_20w_std71

A = Fitted(x1, x2, y1, y2, sd1, sd2, [0, 0, 0], [1, 6, 6])
A.fit_one()


time_all0 = Avrami_data.T_time_all0
survival_all0 = Avrami_data.T_survival_all0
std_all0 = Avrami_data.T_std_all0

time_all1 = Avrami_data.T_time_all1
survival_all1 = Avrami_data.T_survival_all1
std_all1 = Avrami_data.T_std_all1

# """
for i in range(len(time_all0)):
    x1 = [x * 7/365 for x in time_all0[i]]
    x2 = [x * 7/365 for x in time_all1[i]]

    y1 = [1-x for x in survival_all0[i]]
    y2 = [1-x for x in survival_all1[i]]

    sd1 = std_all0[i]
    sd2 = std_all1[i]

    A = Fitted(x1, x2, y1, y2, sd1, sd2)
    A.final_fit_sd()
    A.plot_data_sd()
# """

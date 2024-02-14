import Avrami_data
from Avrami_fits_class import Fitted

x1 = Avrami_data.N_non_mut_x_years
x2 = Avrami_data.N_mut_x_years

y1 = Avrami_data.N_non_mut_canc
y2 = Avrami_data.N_mut_canc

sd1 = Avrami_data.N_non_mut_sd
sd2 = Avrami_data.N_mut_sd

A = Fitted(x1, x2, y1, y2, sd1, sd2)
A.final_fit_sd()
A.plot_data_sd()
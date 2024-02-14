from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


class Fitted(object):
    def __init__(self, x1, x2, y1, y2, sd1=[], sd2=[]):

        self.x1 = x1
        self.y1 = y1
        self.sd1 = sd1

        self.x2 = x2
        self.y2 = y2
        self.sd2 = sd2

        self.comboX = np.append(self.x1, self.x2)
        self.comboY = np.append(self.y1, self.y2)
        self.comboSD = np.append(self.sd1, self.sd2)

        self.popt1, self.pcov1, self.perr1, self.r_squared1 = [], [], [], 0
        self.popt2, self.pcov2, self.perr2, self.r_squared2 = [], [], [], 0
        self.popt3, self.pcov3, self.perr3, self.r_squared31, self.r_squared32 = [], [], [], 0, 0

    def function(self, x, a, n):
        return 1 - np.exp(-a * np.power(x, n))

    def fun1(self, x, a, n1, n2):
        return 1 - np.exp(-a * np.power(x, n1))

    def fun2(self, x, a, n1, n2):
        return 1 - np.exp(-a * np.power(x, n2))

    def comboFunc(self, comboData, a, n1, n2):
        extract1 = comboData[:len(self.y1)]  # first data
        extract2 = comboData[len(self.y1):]  # second data
        result1 = self.fun1(extract1, a, n1, n2)
        result2 = self.fun2(extract2, a, n1, n2)
        return np.append(result1, result2)

    def final_fit_sd(self, bound_low=[0, 0, 0], bound_high=[2, 8, 8]):

        self.popt1, self.pcov1 = curve_fit(self.function, self.x1, self.y1, sigma=self.sd1,
                                 bounds=([bound_low[0], bound_low[1]], [bound_high[0], bound_high[1]]))
        self.perr1 = np.sqrt(np.diag(self.pcov1))  # standard deviation error

        residuals = self.y1 - self.function(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        self.r_squared1 = 1 - (ss_res / ss_tot)

        self.popt2, self.pcov2 = curve_fit(self.function, self.x2, self.y2, sigma=self.sd2,
                                           bounds=([bound_low[0], bound_low[2]], [bound_high[0], bound_high[2]]))
        self.perr2 = np.sqrt(np.diag(self.pcov2))  # standard deviation error

        residuals = self.y2 - self.function(self.x2, *self.popt2)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared2 = 1 - (ss_res / ss_tot)

        self.popt3, self.pcov3 = curve_fit(self.comboFunc, self.comboX, self.comboY, sigma=self.comboSD, bounds=(bound_low, bound_high))
        self.perr3 = np.sqrt(np.diag(self.pcov3))  # standard deviation error - nie zawsze prawdziwe

        print(self.popt3)

        residuals = self.y1 - self.fun1(self.x1, *self.popt3)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        self.r_squared31 = 1 - (ss_res / ss_tot)

        residuals = self.y2 - self.fun2(self.x2, *self.popt3)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared32 = 1 - (ss_res / ss_tot)

    def plot_data_sd(self):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        ax[0].plot(self.x1, [self.function(x, self.popt1[0], self.popt1[1]) for x in self.x1], 'g--',
                   label='fit: a=%.5f, n=%.3f' % (self.popt1[0], self.popt1[1]))
        ax[0].errorbar(self.x1, self.y1, yerr=self.sd1, fmt='o', label="data")
        ax[0].legend()
        ax[0].set_xlabel("time (years)")
        ax[0].set_ylabel("carcinoma probability genotype 0 - BRCA1 +/+")

        ax[0].plot(self.x2, [self.function(x, self.popt2[0], self.popt2[1]) for x in self.x2], 'g--',
                   label='fit: a=%.5f, n=%.3f' % (self.popt2[0], self.popt2[1]))
        ax[0].errorbar(self.x2, self.y2, yerr=self.sd2, fmt='o', label="data")
        ax[0].legend()
        ax[0].set_xlabel("time (years)")
        ax[0].set_ylabel("carcinoma probability genotype 1 - BRCA1 +/+")

        plt.show()

    def return_data(self):
        print(self.popt1, self.pcov1, self.perr1, self.r_squared1)

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np


class Fitted(object):
    def __init__(self, x1, y1, sd1=[], bound_low=[0, 0], bound_high=[2, 10]):
        self.x1 = x1
        self.y1 = y1
        self.y2 = [np.log(-np.log(1 - y)) for y in self.y1]

        self.sd1 = sd1
        self.sd2 = [np.abs(sd1[i]/(np.log(1-y1[i])*(1 - y1[i]))) for i in range(len(sd1))]

        self.bound_low = bound_low
        self.bound_high = bound_high

        self.popt1, self.pcov1, self.perr1, self.r_squared1 = [], [], [], 0
        self.popt_ln, self.pcov_ln, self.perr_ln, self.r_squared_ln = [], [], [], 0

    def function(self, x, a, n):
        return 1 - np.exp(-a * np.power(x, n))

    def function_ln(self, x, a, n):
        return np.log(a) + n * np.log(x)

    def fit_both(self):
        self.popt1, self.pcov1 = curve_fit(self.function, self.x1, self.y1, bounds=(self.bound_low, self.bound_high))
        self.perr1 = np.sqrt(np.diag(self.pcov1))

        residuals = self.y1 - self.function(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        self.r_squared1 = 1 - (ss_res / ss_tot)

        residuals = self.y2 - self.function_ln(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        print(1 - (ss_res / ss_tot))

        self.popt_ln, self.pcov_ln = curve_fit(self.function_ln, self.x1, self.y2, bounds=(self.bound_low, self.bound_high))
        self.perr_ln = np.sqrt(np.diag(self.pcov_ln))  # standard deviation error

        residuals = self.y2 - self.function_ln(self.x1, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared_ln = 1 - (ss_res / ss_tot)

        residuals = self.y1 - self.function(self.x1, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        print(1 - (ss_res / ss_tot))

    def fit_both_sd(self):
        self.popt1, self.pcov1 = curve_fit(self.function, self.x1, self.y1, sigma=self.sd1, bounds=(self.bound_low, self.bound_high))
        self.perr1 = np.sqrt(np.diag(self.pcov1))

        residuals = self.y1 - self.function(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        self.r_squared1 = 1 - (ss_res / ss_tot)

        residuals = self.y2 - self.function_ln(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        print(1 - (ss_res / ss_tot))

        self.popt_ln, self.pcov_ln = curve_fit(self.function_ln, self.x1, self.y2, sigma=self.sd2, bounds=(self.bound_low, self.bound_high))
        self.perr_ln = np.sqrt(np.diag(self.pcov_ln))  # standard deviation error

        residuals = self.y2 - self.function_ln(self.x1, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared_ln = 1 - (ss_res / ss_tot)

        residuals = self.y1 - self.function(self.x1, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        print(1 - (ss_res / ss_tot))

    def plot_data(self):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        ax[0].scatter(self.x1, self.y1, label="data1", color='orange')
        ax[0].plot(self.x1, [self.function(x, self.popt1[0], self.popt1[1]) for x in self.x1], '--',
                   label='fit norm: a=%.5f, n=%.3f' % (self.popt1[0], self.popt1[1]), color='orange')

        ax[0].plot(self.x1, [self.function(x, self.popt_ln[0], self.popt_ln[1]) for x in self.x1], '--',
                   label='fit ln: a=%.5f, n=%.3f' % (self.popt_ln[0], self.popt_ln[1]), color='blue')
        ax[0].legend()
        ax[0].set_ylabel("carcinoma probability")
        ax[0].set_xlabel("time (years)")

        ax[1].scatter(self.x1, self.y2, label="data1", color='blue')
        ax[1].plot(self.x1, [self.function_ln(x, self.popt_ln[0], self.popt_ln[1]) for x in self.x1], '--',
                   label='fit ln: a=%.5f, n=%.3f' % (self.popt_ln[0], self.popt_ln[1]), color='blue')

        ax[1].plot(self.x1, [self.function_ln(x, self.popt1[0], self.popt1[1]) for x in self.x1], '--',
                   label='fit norm: a=%.5f, n=%.3f' % (self.popt1[0], self.popt1[1]), color='orange')

        ax[1].set_xscale('log')
        ax[1].legend()
        ax[1].set_xlabel("time (years)")
        ax[1].set_ylabel("zlogarytmowane")

        plt.show()

    def plot_data_sd(self):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))

        ax[0].scatter(self.x1, self.y1, label="data1", color='orange')
        ax[0].plot(self.x1, [self.function(x, self.popt1[0], self.popt1[1]) for x in self.x1], '--',
                   label='fit norm: a=%.5f, n=%.3f' % (self.popt1[0], self.popt1[1]), color='orange')

        ax[0].plot(self.x1, [self.function(x, self.popt_ln[0], self.popt_ln[1]) for x in self.x1], '--',
                   label='fit ln: a=%.5f, n=%.3f' % (self.popt_ln[0], self.popt_ln[1]), color='blue')
        ax[0].errorbar(self.x1, self.y1, yerr=self.sd1, fmt='o', ecolor='orange', color='orange')
        ax[0].legend()
        ax[0].set_ylabel("carcinoma probability")
        ax[0].set_xlabel("time (years)")

        ax[1].scatter(self.x1, self.y2, label="data1", color='blue')
        ax[1].plot(self.x1, [self.function_ln(x, self.popt_ln[0], self.popt_ln[1]) for x in self.x1], '--',
                   label='fit ln: a=%.5f, n=%.3f' % (self.popt_ln[0], self.popt_ln[1]), color='blue')

        ax[1].plot(self.x1, [self.function_ln(x, self.popt1[0], self.popt1[1]) for x in self.x1], '--',
                   label='fit norm: a=%.5f, n=%.3f' % (self.popt1[0], self.popt1[1]), color='orange')

        ax[1].errorbar(self.x1, self.y2, yerr=self.sd2, fmt='o', ecolor='orange', color='orange')
        ax[1].set_xscale('log')
        ax[1].legend()
        ax[1].set_xlabel("time (years)")
        ax[1].set_ylabel("zlogarytmowane")

        plt.show()
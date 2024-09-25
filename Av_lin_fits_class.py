from scipy.optimize import curve_fit
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

from scipy.stats import linregress


class Fitted(object):
    def __init__(self, x1, y1, sd1=[], bound_low=[0, 0], bound_high=[2, 10]):
        self.x1 = x1
        self.x2 = [np.log(x) for x in self.x1]
        self.y1 = y1
        self.y2 = [np.log(-np.log(1 - y)) for y in self.y1]

        self.sd1 = sd1
        self.sd2 = [np.abs(sd1[i]/(np.log(1-y1[i])*(1 - y1[i]))) for i in range(len(sd1))]

        self.bound_low = bound_low
        self.bound_high = bound_high

        self.popt1, self.pcov1, self.perr1, self.r_squared1 = [], [], [], 0
        self.popt_ln, self.pcov_ln, self.perr_ln, self.r_squared_ln = [], [], [], 0
        self.popt_linear, self.pcov_linear, self.perr_linear, self.r_squared_linear = [], [], [], 0

        self.pear = 0

    def function(self, x, a, n):
        return 1 - np.exp(-a * np.power(x, n))

    def function_ln(self, x, a, n):
        return np.log(a) + n * np.log(x)

    def function_linear(self, x, a, n):
        return a + n * np.array(x)

    def fit_both(self):

        self.pear = pearsonr(self.x2, self.y2)
        print(self.pear)

        # dopasowanie do funkcji Avramiego 1 - np.exp(-a * np.power(x, n))
        self.popt1, self.pcov1 = curve_fit(self.function, self.x1, self.y1, bounds=(self.bound_low, self.bound_high))
        self.perr1 = np.sqrt(np.diag(self.pcov1))
        print(self.popt1)

        residuals = self.y1 - self.function(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        self.r_squared1 = 1 - (ss_res / ss_tot)

        # dopasowanie trochę liniowe np.log(a) + n * np.log(x)
        self.popt_ln, self.pcov_ln = curve_fit(self.function_ln, self.x1, self.y2, bounds=(self.bound_low, self.bound_high))
        self.perr_ln = np.sqrt(np.diag(self.pcov_ln))  # standard deviation error
        print(self.popt_ln)

        residuals = self.y2 - self.function_ln(self.x1, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared_ln = 1 - (ss_res / ss_tot)

        # dopasowanie całkiem liniowe a + n * x
        self.popt_linear, self.pcov_linear = curve_fit(self.function_linear, self.x2, self.y2)
        self.perr_linear = np.sqrt(np.diag(self.pcov_linear))  # standard deviation error
        print(self.popt_linear, np.exp(self.popt_linear[0]))

        residuals = self.y2 - self.function_linear(self.x2, *self.popt_linear)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared_linear = 1 - (ss_res / ss_tot)

        print(stats.linregress([np.log(x) for x in self.x1], self.y2))
        print("\n\n")

    def fit_both_sd(self):
        self.pear = pearsonr(self.x2, self.y2)
        print(self.pear)

        self.popt1, self.pcov1 = curve_fit(self.function, self.x1, self.y1, sigma=self.sd1, bounds=(self.bound_low, self.bound_high))
        self.perr1 = np.sqrt(np.diag(self.pcov1))
        print("1 - np.exp(-a * np.power(x, n))\t", self.popt1, self.perr1)

        residuals = self.y1 - self.function(self.x1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y1 - np.mean(self.y1)) ** 2)
        self.r_squared1 = 1 - (ss_res / ss_tot)

        """
        print(self.y1)
        print([self.function(x, *self.popt1) for x in self.x1])
        print(np.sum(self.y1) / np.sum([self.function(x, *self.popt1) for x in self.x1]) * np.array([self.function(x, *self.popt1) for x in self.x1]))

        chi_square_test_statistic, p_value = stats.chisquare(self.y1, np.sum(self.y1) / np.sum([self.function(x, *self.popt1) for x in self.x1]) * np.array([self.function(x, *self.popt1) for x in self.x1]))
        r_squared1 = r2_score(self.y1, self.function(self.x1, *self.popt1))
        print('chi_square_test_statistic is : ' + str(chi_square_test_statistic))
        print('p_value : ' + str(p_value))
        print(stats.chi2.ppf(1 - 0.05, df=len(self.y1) - 1)) 
        print('r^2 : ' + str(r_squared1))"""

        self.popt_ln, self.pcov_ln = curve_fit(self.function_ln, self.x1, self.y2, sigma=self.sd2, bounds=(self.bound_low, self.bound_high))
        self.perr_ln = np.sqrt(np.diag(self.pcov_ln))  # standard deviation error
        print("np.log(a) + n * np.log(x)\t\t", self.popt_ln, self.perr_ln)

        residuals = self.y2 - self.function_ln(self.x1, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared_ln = 1 - (ss_res / ss_tot)

        self.popt_linear, self.pcov_linear = curve_fit(self.function_linear, self.x2, self.y2, sigma=self.sd2)
        self.perr_linear = np.sqrt(np.diag(self.pcov_linear))  # standard deviation error
        print("a + n * np.array(x)\t\t\t\t", self.popt_linear, np.exp(self.popt_linear[0]), self.perr_linear, np.exp(self.popt_linear[0])*self.perr_linear[0])

        residuals = self.y2 - self.function_linear(self.x2, *self.popt_linear)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y2 - np.mean(self.y2)) ** 2)
        self.r_squared_linear = 1 - (ss_res / ss_tot)

        print(stats.linregress(self.x2, self.y2))

        """chi_square_test_statistic, p_value = stats.chisquare(self.y1, np.sum(self.y1) / np.sum(
            [self.function(x, *self.popt_ln) for x in self.x1]) * np.array(
            [self.function(x, *self.popt_ln) for x in self.x1]))
        r_squared1 = r2_score(self.y1, self.function(self.x1, *self.popt_ln))
        print('chi_square_test_statistic is : ' + str(chi_square_test_statistic))
        print('p_value : ' + str(p_value))
        print(stats.chi2.ppf(1 - 0.05, df=len(self.y1) - 1))
        print('r^2 : ' + str(r_squared1))"""

        print("\n\n")

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

    def plot_data_sd(self, title=""):
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(title)
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

        ax[1].errorbar(self.x1, self.y2, yerr=self.sd2, fmt='o', ecolor='blue', color='blue')
        ax[1].set_xscale('log')
        ax[1].legend()
        ax[1].set_xlabel("time (years)")
        ax[1].set_ylabel("zlogarytmowane")

        plt.show()

    def toasty(self):
        dane1 = {'param': self.popt1, 'cov': self.pcov1, 'err': self.perr1, 'r2': self.r_squared1}
        dane2 = {'param': self.popt_ln, 'cov': self.pcov_ln, 'err': self.perr_ln, 'r2': self.r_squared_ln}
        return dane1, dane2


class Fitted_one_log(object):
    def __init__(self, x, y, sd=[]):
        self.x = [np.log(x1) for x1 in x]
        self.y = [np.log(-np.log(1 - y1)) for y1 in y]
        self.sd = [np.abs(sd[i]/(np.log(1-y[i])*(1 - y[i]))) for i in range(len(sd))]

        self.popt_ln, self.pcov_ln, self.perr_ln, self.r_squared_ln = [], [], [], 0

    def function_ln(self, x, a, n):
        return a + n * np.array(x)

    def fit_both_sd(self):
        self.popt_ln, self.pcov_ln = curve_fit(self.function_ln, self.x, self.y, sigma=self.sd)
        self.perr_ln = np.sqrt(np.diag(self.pcov_ln))  # standard deviation error

        self.pear = pearsonr(self.x, self.y)
        print(self.pear)
        print(stats.linregress(self.x, self.y))

        residuals = self.y - self.function_ln(self.x, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y - np.mean(self.y)) ** 2)
        self.r_squared_ln = 1 - (ss_res / ss_tot)
        return self.popt_ln, self.pcov_ln, self.perr_ln, self.r_squared_ln

    def fit_both(self):

        self.popt_ln, self.pcov_ln = curve_fit(self.function_ln, self.x, self.y)
        self.perr_ln = np.sqrt(np.diag(self.pcov_ln))  # standard deviation error

        self.pear = pearsonr(self.x, self.y)
        # print(self.pear)

        residuals = self.y - self.function_ln(self.x, *self.popt_ln)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.y - np.mean(self.y)) ** 2)
        self.r_squared_ln = 1 - (ss_res / ss_tot)
        return self.popt_ln, self.pcov_ln, self.perr_ln, self.r_squared_ln

    def toasty(self):
        dane1 = {'param': self.popt_ln, 'cov': self.pcov_ln, 'err': self.perr_ln, 'r2': self.r_squared_ln}
        return dane1



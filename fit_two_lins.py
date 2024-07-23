from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import ks_2samp
from scipy.stats import linregress


class Fitted(object):
    def __init__(self, x1, y1, sd1=[], window=10, start=0, stop=0):
        self.xx = x1
        self.yy = y1

        self.sd1 = sd1

        self.popt1, self.pcov1, self.perr1, self.r_squared1 = [], [], [], 0
        self.popt2, self.pcov2, self.perr2, self.r_squared2 = [], [], [], 0

        self.xx1 = []
        self.yy1 = []

        self.xx2 = []
        self.yy2 = []

    def function(self, x, a, n):
        return a + n * np.array(x)

    def fit_both(self):

        # a = n, b = ln(k) => k = e^b

        """
        self.popt1, self.pcov1 = curve_fit(self.function, self.xx1, self.yy1)
        self.perr1 = np.sqrt(np.diag(self.pcov1))
        print("n = " + str(self.popt1[1]) + ",\t k = " + str(np.exp(self.popt1[0])) + ",\t b = " + str(self.popt1[0]))
        print("n = " + str(self.perr1[1]) + ",\t b = " + str(self.perr1[0]))

        residuals = self.yy1 - self.function(self.xx1, *self.popt1)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.yy1 - np.mean(self.yy1)) ** 2)
        self.r_squared1 = 1 - (ss_res / ss_tot)"""

        do = stats.linregress(self.xx1, self.yy1)
        a1, b1, r2_1, pval1, ua1, ub1 = do.slope, do.intercept, do.rvalue, do.pvalue, do.stderr, do.intercept_stderr

        # dopasowanie trochę liniowe np.log(a) + n * np.log(x)
        """
        self.popt2, self.pcov2 = curve_fit(self.function, self.xx2, self.yy2)
        self.perr2 = np.sqrt(np.diag(self.pcov2))  # standard deviation error
        print("n = " + str(self.popt2[1]) + ",\t k = " + str(np.exp(self.popt2[0])) + ",\t b = " + str(self.popt2[0]))
        print("n = " + str(self.perr2[1]) + ",\t b = " + str(self.perr2[0]))

        residuals = self.yy2 - self.function(self.xx2, *self.popt2)
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((self.yy2 - np.mean(self.yy2)) ** 2)
        self.r_squared2 = 1 - (ss_res / ss_tot)
        print(self.r_squared2)"""

        print(stats.linregress(self.xx2, self.yy2))
        do = stats.linregress(self.xx2, self.yy2)
        a2, b2, r2_2, pval2, ua2, ub2 = do.slope, do.intercept, do.rvalue, do.pvalue, do.stderr, do.intercept_stderr
        print("\n\n")
        return a1, b1, r2_1, pval1, ua1, ub1, a2, b2, r2_2, pval2, ua2, ub2

    def fitt(self, start=20, okno=10):

        self.xx1 = self.xx[:start]
        self.yy1 = self.yy[:start]

        self.xx2 = self.xx[start+okno-1:]
        self.yy2 = self.yy[start+okno-1:]

        return self.fit_both()


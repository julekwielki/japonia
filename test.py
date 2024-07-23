
from Avrami_fits_class import Fitted_lin
from Avrami_fits_class import Fitted_one
from fit_two_lins import Fitted
import matplotlib.pyplot as plt
import random
import pandas as pd
import scipy.stats as stats
import numpy as np

cols = ['C:E', 'F:H', 'I:K', 'L:N', 'O:Q', 'R:T', 'U:W', 'X:Z', 'AA:AC', 'AD:AF', 'AG:AI', 'AJ:AL', 'AM:AO', 'AP:AR', 'AS:AU', 'AV:AX', 'AY:BA']

for c in cols:
    df = pd.read_excel('dane_cancer_mutation_age.xlsx', usecols=c)  # , nrows=10
    xx = []
    yy = []
    print(df.columns[0])
    print(df.shape[0])
    for x in range(df.shape[0]-1):
        if df.iloc[x+1, 1] > 0:
            xx.append(df.iloc[x+1, 1]/100)
            yy.append(df.iloc[x+1, 2])
    xx.pop(0)
    yy.pop(0)
    yy = [x for _, x in sorted(zip(xx, yy))]
    xx = sorted(xx)

    unique_hist, counts_hist = np.unique(xx, return_counts=True)

    xx2 = [sum(counts_hist[:x]) for x in range(len(counts_hist))]
    plt.scatter(unique_hist, counts_hist)
    plt.show()
    # plt.xlabel("wiek")
    # plt.ylabel("liczba pacjentów")
    # plt.show()

    # fits = Fitted_one(unique_hist, xx2)
    # fits.final_fit()
    # fits.plot_data(df.columns[0])

cols = ['C:E', 'F:H', 'I:K', 'L:N', 'O:Q']

for c in cols:
    df = pd.read_excel('dane_cancer_mutation_age.xlsx', sheet_name='(MSK, Nat Genet 2020)', usecols=c)  # , nrows=10
    xx = []
    yy = []
    print(df.columns[0])
    print(df.shape[0])
    for x in range(df.shape[0]-1):
        if df.iloc[x+1, 1] > 0:
            xx.append(df.iloc[x+1, 1]/100)
            yy.append(df.iloc[x+1, 2])
    xx.pop(0)
    yy.pop(0)
    yy = [x for _, x in sorted(zip(xx, yy))]

    unique_hist, counts_hist = np.unique(xx, return_counts=True)

    xx2 = [sum(counts_hist[:x]) for x in range(len(counts_hist))]
    plt.scatter(unique_hist, xx2)
    plt.show()

    """fits = Fitted_one(unique_hist, xx2)
    fits.final_fit()
    fits.plot_data(df.columns[0])"""


"""
cols = ['C:E', 'F:H', 'I:K', 'L:N', 'O:Q', 'R:T', 'U:W', 'X:Z', 'AA:AC', 'AD:AF', 'AG:AI', 'AJ:AL', 'AM:AO', 'AP:AR', 'AS:AU', 'AV:AX', 'AY:BA']

for c in cols:
    df = pd.read_excel('dane_cancer_mutation_age.xlsx', usecols=c)  # , nrows=10
    xx = []
    yy = []
    print(df.columns[0])
    print(df.shape[0])
    for x in range(df.shape[0]-1):
        if df.iloc[x+1, 1] > 0:
            xx.append(df.iloc[x+1, 1])
            yy.append(df.iloc[x+1, 2])
    xx.pop(0)
    yy.pop(0)
    yy = [x for _, x in sorted(zip(xx, yy))]
    xx = sorted(xx)
    print(max(yy))
    plt.scatter(xx, yy)
    # plt.show()
    do = stats.linregress(xx, yy)
    a1, b1, r2_1, pval1, ua1, ub1 = do.slope, do.intercept, do.rvalue, do.pvalue, do.stderr, do.intercept_stderr

    fits = Fitted_lin(xx, yy)
    fits.final_fit()
    fits.plot_data(df.columns[0])
# """

"""
cols = ['C:E', 'F:H', 'I:K', 'L:N', 'O:Q']

for c in cols:
    df = pd.read_excel('dane_cancer_mutation_age.xlsx', sheet_name='(MSK, Nat Genet 2020)', usecols=c)  # , nrows=10
    xx = []
    yy = []
    print(df.columns[0])
    print(df.shape[0])
    for x in range(df.shape[0]-1):
        if df.iloc[x+1, 1] > 0:
            xx.append(df.iloc[x+1, 1])
            yy.append(df.iloc[x+1, 2])
    xx.pop(0)
    yy.pop(0)
    yy = [x for _, x in sorted(zip(xx, yy))]
    xx = sorted(xx)
    print(max(yy))
    plt.scatter(xx, yy)
    # plt.show()
    do = stats.linregress(xx, yy)
    a1, b1, r2_1, pval1, ua1, ub1 = do.slope, do.intercept, do.rvalue, do.pvalue, do.stderr, do.intercept_stderr

    fits = Fitted_lin(xx, yy)
    fits.final_fit()
    fits.plot_data(df.columns[0])
#"""

"""
plt.plot(xx, yy)
plt.show()
start = 20
okno = 10
toasty = Fitted(xx, yy, 20, 10)

a1, b1, r2_1, pval1, ua1, ub1, a2, b2, r2_2, pval2, ua2, ub2 = toasty.fitt()

x_cross = (b2 - b1)/(a1-a2)
print(x_cross)

xx_fit1 = []
xx_fit2 = []
yy_fit1 = []
yy_fit2 = []
for x in xx:
    if x <= start:
        xx_fit1.append(x)
        yy_fit1.append(a1 * x + b1)
    elif start < x < start + okno:
        xx_fit1.append(x)
        yy_fit1.append(a1 * x + b1)

        xx_fit2.append(x)
        yy_fit2.append(a2 * x + b2)
    else:
        xx_fit2.append(x)
        yy_fit2.append(a2 * x + b2)

plt.plot(xx, yy)
plt.plot(xx_fit1, yy_fit1)
plt.plot(xx_fit2, yy_fit2)
plt.show()"""


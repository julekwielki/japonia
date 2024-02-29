import numpy as np
import random
import matplotlib.pyplot as plt
from Av_lin_fits_class import Fitted

length = 200
p = 0.001
ziarna = [[], []]
ziarna_t = []
data = np.zeros((length, length))
av = []
time = []


def add_ziarna(t):
    for i in range(length):
        for j in range(length):
            if data[i][j] == 0:
                if random.random() < p:
                    data[i][j] = 1
                    ziarna[0].append(i)
                    ziarna[1].append(j)
                    ziarna_t.append(t)


def ziarno_xy(x, y, t):
    data[x][y] = 1
    ziarna[0].append(x)
    ziarna[1].append(y)
    ziarna_t.append(t)


# add_ziarna(0)
# print(len(ziarna_t))
# time.append(0)
# av.append(len(ziarna_t))
ziarno_xy(0, 0, 0)
# """
if len(ziarna_t) > 0:
    for t in range(1, 50):
        print(t)
        for i in range(length):
            for j in range(length):
                if data[i][j] == 0:
                    for x in range(len(ziarna[0])):
                        dist = (ziarna[0][x] - i)**2 + (ziarna[1][x] - j)**2
                        if dist <= (ziarna_t[x] - t)**2:
                            data[i][j] = 1
        time.append(t)
        a = 0
        for i in range(length):
            for j in range(length):
                if data[i][j] != 0:
                    a += 1
        av.append(a)

    avr_log, time_log, time2, av2 = [], [], [], []
    for i in range(len(av)):
        if 0 < av[i] < length**2:
            avr_log.append(np.log(-np.log(1 - av[i] / length**2)))
            time2.append(time[i]+1 / 10)
            time_log.append(np.log(time[i]+1))
            av2.append(av[i] / length**2)
    # fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    # ax[0].scatter(time, av)
    # ax[1].scatter(time_log, avr_log)
    #plt.show()

    A = Fitted(time2, av2)
    A.fit_both()
    A.plot_data()
# """  # brak nowych ziaren ze wzrostem ziaren

"""
for t in range(1, 50):
    print(t)

    add_ziarna(t)
    for i in range(length):
        for j in range(length):
            if data[i][j] == 0:
                for x in range(len(ziarna[0])):
                    dist = (ziarna[0][x] - i)**2 + (ziarna[1][x] - j)**2
                    if dist <= (ziarna_t[x] - t)**2:
                        data[i][j] = 1
    time.append(t)
    a = 0
    for i in range(length):
        for j in range(length):
            if data[i][j] != 0:
                a += 1
    av.append(a)

avr_log, time_log, time2, av2 = [], [], [], []
for i in range(len(av)):
    if 0 < av[i] < length**2:
        avr_log.append(np.log(-np.log(1 - av[i] / length**2)))
        time2.append(time[i]+1 / 50)
        time_log.append(np.log(time[i]+1))
        av2.append(av[i] / length**2)
# fig, ax = plt.subplots(1, 2, figsize=(12, 5))
# ax[0].scatter(time, av)
# ax[1].scatter(time_log, avr_log)
# plt.show()

A = Fitted(time2, av2)
A.fit_both()
A.plot_data()
# """  # nowe ziarna ze wzrostem ziaren


"""
for t in range(1, 50):
    print(t)
    add_ziarna(t)
    time.append(t)
    a = 0
    for i in range(length):
        for j in range(length):
            if data[i][j] != 0:
                a += 1
    av.append(a)

avr_log, time_log, time2, av2 = [], [], [], []
for i in range(len(av)):
    if 0 < av[i] < length**2:
        avr_log.append(np.log(-np.log(1 - av[i] / length**2)))
        time2.append(time[i]+1 / 50)
        time_log.append(np.log(time[i]+1))
        av2.append(av[i] / length**2)
# fig, ax = plt.subplots(1, 2, figsize=(12, 5))
# ax[0].scatter(time, av)
# ax[1].scatter(time_log, avr_log)
# plt.show()

A = Fitted(time2, av2)
A.fit_both()
A.plot_data()
# """  # nowe ziarna bez wzrostu ziaren


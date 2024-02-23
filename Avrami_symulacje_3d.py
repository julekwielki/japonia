import numpy as np
import random
import matplotlib.pyplot as plt
from Av_lin_fits_class import Fitted

length = 80
p = 0.00005
ziarna = [[], [], []]
ziarna_t = []
data = np.zeros((length, length, length))
av = []
time = []


def add_ziarna(t):
    for i in range(length):
        for j in range(length):
            for k in range(length):
                if data[i][j][k] == 0:
                    if random.random() < p:
                        data[i][j][k] = 1
                        ziarna[0].append(i)
                        ziarna[1].append(j)
                        ziarna[2].append(k)
                        ziarna_t.append(t)


add_ziarna(0)
print(len(ziarna_t))
# time.append(0)
# av.append(len(ziarna_t))
"""
if len(ziarna_t) > 0:
    for t in range(1, 100):
        print(t)

        for i in range(length):
            for j in range(length):
                for k in range(length):
                    if data[i][j][k] == 0:
                        for x in range(len(ziarna[0])):
                            dist = (ziarna[0][x] - i) ** 2 + (ziarna[1][x] - j) ** 2 + (ziarna[2][x] - k) ** 2
                            if dist <= (ziarna_t[x] - t) ** 2:
                                data[i][j][k] = 1

        time.append(t)
        a = 0
        for i in range(length):
            for j in range(length):
                for k in range(length):
                    if data[i][j][k] != 0:
                        a += 1
        av.append(a)
        if a == length**3:
            break

    avr_log, time_log, time2, av2 = [], [], [], []
    for i in range(len(av)):
        if 0 < av[i] < length ** 3:
            avr_log.append(np.log(-np.log(1 - av[i] / length ** 3)))
            time2.append(time[i] + 1 / 50)
            time_log.append(np.log(time[i] + 1))
            av2.append(av[i] / length ** 3)
    # fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    # ax[0].scatter(time, av)
    # ax[1].scatter(time_log, avr_log)
    # plt.show()

    A = Fitted(time2, av2)
    A.fit_both()
    A.plot_data()
# """  # brak nowych ziaren ze wzrostem ziaren

"""
for t in range(1, 100):
    print(t)

    add_ziarna(t)
    for i in range(length):
        for j in range(length):
            for k in range(length):
                if data[i][j][k] == 0:
                    for x in range(len(ziarna[0])):
                        dist = (ziarna[0][x] - i) ** 2 + (ziarna[1][x] - j) ** 2 + (ziarna[2][x] - k) ** 2
                        if dist <= (ziarna_t[x] - t) ** 2:
                            data[i][j][k] = 1

    time.append(t)
    a = 0
    for i in range(length):
        for j in range(length):
            for k in range(length):
                if data[i][j][k] != 0:
                    a += 1
    av.append(a)
    if a == length**3:
        print("breaky")
        break

avr_log, time_log, time2, av2 = [], [], [], []
for i in range(len(av)):
    if 0 < av[i] < length**3:
        avr_log.append(np.log(-np.log(1 - av[i] / length**3)))
        time2.append(time[i]+1 / 50)
        time_log.append(np.log(time[i]+1))
        av2.append(av[i] / length**3)
# fig, ax = plt.subplots(1, 2, figsize=(12, 5))
# ax[0].scatter(time, av)
# ax[1].scatter(time_log, avr_log)
# plt.show()

A = Fitted(time2, av2)
A.fit_both()
A.plot_data()
# """  # nowe ziarna ze wzrostem ziaren


# """
for t in range(1, 50):
    print(t)
    add_ziarna(t)
    time.append(t)
    a = 0
    for i in range(length):
        for j in range(length):
            for k in range(length):
                if data[i][j][k] != 0:
                    a += 1
    av.append(a)
    if a == length**3:
        break

avr_log, time_log, time2, av2 = [], [], [], []
for i in range(len(av)):
    if 0 < av[i] < length**3:
        avr_log.append(np.log(-np.log(1 - av[i] / length**3)))
        time2.append(time[i]+1 / 50)
        time_log.append(np.log(time[i]+1))
        av2.append(av[i] / length**3)
# fig, ax = plt.subplots(1, 2, figsize=(12, 5))
# ax[0].scatter(time, av)
# ax[1].scatter(time_log, avr_log)
# plt.show()

A = Fitted(time2, av2)
A.fit_both()
A.plot_data()
# """  # nowe ziarna bez wzrostu ziaren


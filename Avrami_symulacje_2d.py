import numpy as np
import random
import matplotlib.pyplot as plt
from Av_lin_fits_class import Fitted

length = 10
p = 0.05
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


add_ziarna(0)

# """
if len(ziarna) > 0:
    for t in range(1, 100):
        print(t)
        # add_ziarna(t)
        for i in range(len(data)):
            if data[i] == 0:

                for j in range(len(ziarna)):
                    if np.abs(ziarna[j] - i) <= np.abs(ziarna_t[j] - t):
                        # print(i, j, ziarna[j], ziarna_t[j], np.abs(ziarna[j] - i), np.abs(ziarna_t[j]-t))
                        data[i] = 1

        time.append(t)
        a = 0
        for i in data:
            if i != 0:
                a += 1
        av.append(a)

    av.pop(0)
    time.pop(0)

    avr_log, time_log, time2, av2 = [], [], [], []

    for i in range(len(av)):
        if 0 < av[i] < len(data):
            avr_log.append(np.log(-np.log(1 - av[i] / len(data))))
            time2.append(time[i] / 10)
            time_log.append(np.log(time[i]))
            av2.append(av[i] / len(data))

    # fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    # ax[0].scatter(time, av)
    # ax[1].scatter(time_log, avr_log)
    # plt.show()

    A = Fitted(time2, av2)
    A.fit_both()
    A.plot_data()
# """  # nowe ziarna lub nie ze wzrostem ziaren


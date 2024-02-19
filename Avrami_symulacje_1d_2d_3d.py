import numpy as np
import random
import matplotlib.pyplot as plt

len = 10
p = 0.01
ziarna = []
ziarna_t = []
data = np.zeros(len)


def fun_ziarna():
    for i in range(data.size):
        if random.random() < p:
            data[i] = 1
            ziarna.append(i)


def add_ziarna():
    for i in range(data.size):
        if data[i] == 0:
            if random.random() < p:
                data[i] = 1
                ziarna.append(i)
                ziarna_t.append(t)

avr = []
time = []
for t in range(20):
    time.append(t)

    if t == 0:
        add_ziarna()

print(data)
print(ziarna)
print(ziarna_t)


"""
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].scatter(time, avr)

avr.pop(0)
time.pop(0)

avr2 = [np.log(-np.log(1 - y/len)) for y in avr]
time2 = [np.log(x) for x in time]

ax[1].scatter(time2, avr2)
plt.show()
"""

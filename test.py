
from fit_two_lins import Fitted
import matplotlib.pyplot as plt
import random

xx = []
yy = []

x1 = 20
x2 = 30
xp = 22

a1 = 5
b1 = -10

a3 = 1
b3 = a1 * xp + b1 - a3 * xp

a2 = (a3 * x2 + b3 - a1 * x1 - b1)/(x2 - x1)
b2 = a1 * x1 + b1 - a2 * x1

print(a1, b1)
print(a2, b2)
print(a3, b3)

for x in range(1, 50):
    if x <= x1:
        xx.append(x)
        yy.append(a1 * x + b1 + random.random() - 0.5)
    elif x1 <= x < x2:
        xx.append(x)
        yy.append(a2 * x + b2 + random.random() - 0.5)
    else:
        xx.append(x)
        yy.append(a3 * x + b3 + random.random() - 0.5)

plt.plot(xx, yy)
plt.show()

toasty = Fitted(xx, yy)

toasty.fitt()
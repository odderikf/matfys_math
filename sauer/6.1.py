#%% setup
import numpy as np
from matplotlib import pyplot as plt


def euler(f, y_0, start, stop, step):
    w_i = y_0
    w = [w_i]
    t = np.arange(start, stop+step, step)  # include endpoint, therefore add step
    for t_i in t[:-1]:
        w_i += step*f(t_i, w_i)
        w.append(w_i)
        
    return t, w


def y_a(t):
    return 0.5*t*t + 1


def y_b(t):
    return np.e ** (t*t*t/3)


def dydt_a(t, y):
    return t


def dydt_b(t, y):
    return t*t*y

#%% run
h = 0.1
a, b = 0., 1.
y_0 = 1.
t, w_a = euler(dydt_a, y_0, a, b, h)
_, w_b = euler(dydt_b, y_0, a, b, h)

#%% print a
print("t:       ", '|', "w:       ", '|', "error:  ")
for i in range(len(t)):
    print(f'{t[i]:+6f} | {w_a[i]:+6f} | {y_a(t[i]) - w_a[i]:+6f}')

#%% print b ØVING
print("t:       ", '|', "w:       ", '|', "error:  ")
for i in range(len(t)):
    print(f'{t[i]:+6f} | {w_b[i]:+6f} | {y_b(t[i]) - w_b[i]:+6f}')

#%% plot a
plt.plot(t, w_a, 'r')
plt.scatter(t, w_a, c='r', s=15)
plt.plot(t, y_a(t), 'b')
plt.show()

#%% plot b
plt.plot(t, w_b, 'r')
plt.scatter(t, w_b, c='r', s=15)
plt.plot(t, y_b(t), 'b')
plt.show()
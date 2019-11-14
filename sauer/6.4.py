# 1b

#%% setup
import numpy as np
from matplotlib import pyplot as plt


def midpoint(f, y_0, start, stop, step): # todo do
    w_i = y_0
    w = [w_i]
    t = np.arange(start, stop+step, step)  # include endpoint, therefore add step
    for t_i in t[:-1]:  # skip last, because i'm using t_i to find w_i+1
        f_i = f(t_i, w_i)
        w_i += step*f(t_i+step/2, w_i + f_i*h/2) / 2
        w.append(w_i)
    return t, w


def RK4(f, y_0, start, stop, step):
    w_i = y_0
    w = [w_i]
    t = np.arange(start, stop+step, step)  # include endpoint, therefore add step
    for t_i in t[:-1]:  # skip last, because i'm using t_i to find w_i+1
        s_i_0 = f(t_i, w_i)
        s_i_1 = f(t_i + h/2, w_i + 0.5*s_i_0)
        s_i_2 = f(t_i + h/2, w_i + 0.5*s_i_1)
        s_i_3 = f(t_i + h, w_i + h*s_i_2)
        w_i += step/6 * (s_i_0 + 2*s_i_1 + 2*s_i_2 + s_i_3)
        w.append(w_i)
    return t, w


def y(t):
    return np.e ** (t*t*t/3)


def dydt(t, y):
    return t*t*y


#%% run
h = 0.1
a, b = 0., 1.
y_0 = 1.
t, w_midpoint = midpoint(dydt, y_0, a, b, h)
t, w_RK4 = RK4(dydt, y_0, a, b, h)

#%% print a
print("t:       ", '|', "w_m:     ", '|', "w_r:     ", '|', "error:  ")
for i in range(len(t)):
    print(f'{t[i]:+6f} | {w_midpoint[i]:+6f} | {w_RK4[i]:+6f} | {y(t[i]) - w_midpoint[i]:+6f}')

#%% plot a
plt.plot(t, w_midpoint, 'r')
plt.scatter(t, w_midpoint, c='r', s=15)
plt.plot(t, y(t), 'b')
plt.show()
plt.plot(t, w_RK4, 'r')
plt.scatter(t, w_RK4, c='r', s=15)
plt.plot(t, y(t), 'b')
plt.show()

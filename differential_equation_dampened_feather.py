import RK45
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np
matplotlib.use('TkAgg')


def dydt(t, y):
    # the equation Ay'' + By' + Cy = 0
    # in context of specifically a block with friction on a spring, with y = displacement
    # so ma + µv + kx = mx'' + µx' + kx = 0
    v, x = y
    F = -g*v - k*x
    a = F/m
    return np.array([a, v])


m = 1
k = 1

y0 = [0, 2]
a = 0
b = 30
h = 0.1
ts = []
ws = []
n = 200
gs = np.linspace(0, 3, n)
for g in gs:

    t, w = RK45.RK45(dydt, y0, a, b, h, 1e-8)
    ts.append(t)
    w = np.array(w)
    ws.append(w[:, 1])

g = np.sqrt(4*m*k)
t_crit, w_crit = RK45.RK45(dydt, y0, a, b, h, 1e-6)
w_crit = np.array(w_crit)[:, 1]

rootparts = [1 - 4*m*k / e**2 for e in gs ]
num_neg_roots = len([e for e in rootparts if e < 0])

red_cs = [(1-i, 0.5*i, 0, 0.6) for i in np.linspace(0, 1, num_neg_roots)]
blue_cs = [(0, 0.5*(1-i), i, 0.6) for i in np.linspace(0, 1, n - num_neg_roots)]

for i in range(n):
    t = ts[i]
    w = ws[i]
    rootpart = rootparts[i]
    print(rootpart, i)
    if rootpart < 0:
        c = red_cs[i]
    else:
        c = blue_cs[i - num_neg_roots]
    plt.plot(t, w, c=c, lw=1)

plt.plot(t_crit, w_crit, c=(0, 1, 0), lw=1.2)
plt.grid(True)
red_patch = mpatches.Patch(color=(1,0,0), label='Dampened oscillation, two complex roots')
green_patch = mpatches.Patch(color=(0,1,0), label='Critical dampening, one real root')
blue_patch = mpatches.Patch(color=(0,0,1), label='Overdampening, two real roots')
plt.legend(handles=[red_patch, green_patch, blue_patch])
plt.show()
plt.savefig('tight.png', dpi=300)
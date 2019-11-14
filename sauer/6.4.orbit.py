#%% setupbytt
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
matplotlib.use('TkAgg')


def RK4(f, y_0, start, stop, step):
    w_i = np.array([i for i in y_0])
    w = [w_i]
    t = np.arange(start, stop+step, step)  # include endpoint, therefore add step
    for t_i in t[:-1]:  # skip last, because i'm using t_i to find w_i+1
        s_i_0 = f(t_i, w_i)
        s_i_1 = f(t_i + h/2, w_i + 0.5*s_i_0)
        s_i_2 = f(t_i + h/2, w_i + 0.5*s_i_1)
        s_i_3 = f(t_i + h, w_i + h*s_i_2)
        w_i = np.array([i for i in w_i]) + (step/6) * (s_i_0 + 2*s_i_1 + 2*s_i_2 + s_i_3)
        w.append(w_i)
    return t, w


def dydt(t, y_v):
    x, y, z = y_v
    return np.array([
        -s*x + s*y,
        -x*z + r*x - y,
        x*y-b*z
    ])


#%% run
h = 0.000001
start, stop = 0, 2
y_0 = (5., 5., 5)
s = 10.
r = 28.
b = 8./3.

t, w = RK4(dydt, y_0, start, stop, h)
w = np.array(w)[::10]

#%% plot

fig = plt.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-15, 15), ylim=(+0, 30))

line1, = axes.plot([], [], 'g', lw=2)


def init():
    """initialize animation"""
    line1.set_data([], [])
    return line1,


def animate(i):
    """perform animation step"""
    line1.set_data(w[0:i, 0], w[0:i, 2])  # xz
    return line1,


anim=animation.FuncAnimation(fig,  # figure to plot in
                             animate,  # function that is called on each frame
                             frames=len(w),  # total number of frames
                             repeat=False,
                             interval=1,
                             blit=True,
                             init_func=init  # initialization
                             )

plt.show()




#%% setup
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import animation
matplotlib.use('TkAgg')


def trapezoid(f, y_0, start, stop, step):
    w_i = np.array([i for i in y_0])
    w = [w_i]
    t = np.arange(start, stop+step, step)  # include endpoint, therefore add step
    for t_i in t[:-1]:  # skip last, because i'm using t_i to find w_i+1
        f_i = f(t_i, w_i)
        w_i = np.array([i for i in w_i]) + step*(f_i + f(t_i+step, w_i + h*f_i)) / 2
        w.append(w_i)
    return t, w


def dydt(t, y):
    x0, y0, vx0, vy0, x1, y1, vx1, vy1, x2, y2, vx2, vy2 = y
    r01sq = (x1 - x0) ** 2 + (y1 - y0) ** 2
    r02sq = (x2 - x0) ** 2 + (y2 - y0) ** 2
    r12sq = (x2 - x1) ** 2 + (y2 - y1) ** 2
    a01 = G / r01sq
    a02 = G / r02sq
    a12 = G / r12sq

    ax0 = m[1]*a01*(x1-x0)/r01sq**0.5 + m[2]*a02*(x2-x0)/r02sq**0.5
    ay0 = m[1]*a01*(y1-y0)/r01sq**0.5 + m[2]*a02*(y2-y0)/r02sq**0.5
    ax1 = m[0]*a01*(x0-x1)/r01sq**0.5 + m[2]*a12*(x2-x1)/r12sq**0.5
    ay1 = m[0]*a01*(y0-y1)/r01sq**0.5 + m[2]*a12*(y2-y1)/r12sq**0.5
    ax2 = m[0]*a02*(x0-x2)/r02sq**0.5 + m[1]*a12*(x1-x2)/r12sq**0.5
    ay2 = m[0]*a02*(y0-y2)/r02sq**0.5 + m[1]*a12*(y1-y2)/r12sq**0.5

    return np.array([
        vx0,
        vy0,
        ax0,
        ay0,

        vx1,
        vy1,
        ax1,
        ay1,

        vx2,
        vy2,
        ax2,
        ay2,
    ])

#%% run
G = 1
h = 0.001
a, b = 0., 180
y_0 = [2., 2., 0.2, -0.2, 0., 0., 0., 0., -2., -2., -0.2, 0.2]  # x, y, vx, vy, three times
m = [0.03, 0.3, 0.03]
t_1, w_1 = trapezoid(dydt, y_0, a, b, h)
G = 1
h = 0.01
a, b = 0., 180
y_0 = [-0.970, 0.243, -0.466, -0.433, 0.970, -0.243, -0.466, -0.433, 0, 0, 2*0.466, 2*0.433]  # x, y, vx, vy, three times
m = [1., 1., 1.]
t_2, w_2 = trapezoid(dydt, y_0, a, b, h)

#%% plot

w = np.array(w_1)

fig = plt.figure()
axes = fig.add_subplot(111, aspect='equal', autoscale_on=False,
                     xlim=(-6, 6), ylim=(-6, 6))

line1, = axes.plot([], [], 'o-g', lw=2)
line2, = axes.plot([], [], 'o-r', lw=2)
line3, = axes.plot([], [], 'o-b', lw=2)


def init():
    """initialize animation"""
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return line1, line2, line3


def animate(i):
    """perform animation step"""
    tail = i-100 if i>100 else 0
    line1.set_data(w[i, 0], w[i, 1])
    line2.set_data(w[i, 4], w[i, 5])
    line3.set_data(w[i, 8], w[i, 9])
    return line1, line2, line3


anim=animation.FuncAnimation(fig,  # figure to plot in
                             animate,  # function that is called on each frame
                             frames=len(w),  # total number of frames
                             repeat=False,
                             interval=1,
                             blit=True,
                             init_func=init  # initialization
                             )

plt.show()
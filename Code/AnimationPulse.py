# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 12:32:41 2018

@author: Lukas
"""
import numpy as np

import matplotlib.pyplot as plt
#from matplotlib.animation import FuncAnimation
from DiffusionEquationBSplines import t,xfine
from matplotlib.animation import FuncAnimation

# =============================================================================
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# =============================================================================

phi = np.loadtxt('phipulse')
L = xfine[-1]


fig, ax = plt.subplots()
ln, = ax.plot([], [], 'r', animated=True,label='Laserpulse')


def init():
    ax.set_xlim(0, L) #max x value of the animation window
    ax.set_ylim(0, 2) #maximum y value the animation can have
    return ln,


def update(frame):
    ln.set_data(xfine,phi[frame])
    return ln,


ani = FuncAnimation(fig, update,init_func=init, blit=True)
plt.xlabel('Depth of Material')
plt.ylabel('Temperature')
plt.legend(loc='best')
plt.show()





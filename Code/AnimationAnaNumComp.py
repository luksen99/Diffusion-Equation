# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 12:32:41 2018
This code is dependent on the Program: DirichletBoundaryConditions
@author: Lukas
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from DiffusionEquationBSplines import t,xfine

phi = np.loadtxt('phi1')
#phi = np.loadtxt('phipulse')
L = xfine[-1]
omega = (np.pi/L)


fig, ax = plt.subplots()
ln, = ax.plot([], [], 'r', animated=True,label='numeric')
ln2, = ax.plot([],[],'k',animated = True,label= 'analytic')

def init():
    ax.set_xlim(0, L) #max x value of the animation window
    ax.set_ylim(0, 2) #maximum y value the animation can have
    return ln, ln2,


def update(frame):
    
    ln2.set_data(xfine,np.sin(omega*xfine)*np.exp(-omega**2*t[frame]))
    ln.set_data(xfine,phi[frame])
    return ln, ln2,


ani = FuncAnimation(fig, update,init_func=init, blit=True)
plt.xlabel('Depth of Material')
plt.ylabel('Temperature')
plt.legend(loc='best')
plt.show()

#


"""
Created on Tue Jun 26 12:32:41 2018
This animation depends on the safed file 'phipulse' from the DirichletBoundaryConditions file
@author: Lukas
"""
import numpy as np
import matplotlib.pyplot as plt
from DiffusionEquationBSplines import t,xfine
from matplotlib.animation import FuncAnimation

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

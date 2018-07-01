
"""
Created on Sat Jun 30 14:52:34 2018
Implementing a function k(phi) = k0 +k1*phi
and solve the heat diffusion equation
@author: Lukas
"""
import numpy as np
import scipy.linalg as spl
from bspline import Bspline
from bspline.splinelab import aptknt

x     = np.linspace(0,10,10)   #space grid to solve for coefficients c
xfine = np.linspace(0,10,1000) #gird to plot
dt    = .05                    #timesteps
t     = np.arange(0,20,dt)     #timegrid 
order = 5                      #order of B-spline

knot_vector = aptknt(x, order)     #kreates knot points with Ghost points
basis = Bspline(knot_vector,order) #Object of basic spline vectors: basis(0) gives 0´th basis spline

#==================Creating 4 Splines: A,A1,B,C ===============================

C = basis.collmat(xfine) #spline to plot/evaluate solution once coefficients c are obtained
C[-1,-1] = 1             #sets the bottom right point to 1 instead of 0

B = basis.collmat(x) #creates the Matrix Bij
B[-1,-1] = 1         #Bottom right point to 1 instead of 0 to avoid singular matrix
Binv = spl.inv(B) 

A  = basis.collmat(x,deriv_order=2) #creates Matrix Aij; 2n order derivative of splines 

A1 = basis.collmat(x,deriv_order=1) #creates Matrix Aij; 2n order derivative of splines 

#================= initial coefficients =======================================
c = np.dot(Binv,np.sin(np.pi/x[-1]*x)) 

#================= Dirichlet Boundary conditions ==============================
bc = np.zeros(np.size(B,0)); bc[0] = 0; bc[-1] = 0

k0 = 1; k1 = 0.2 #coefficints of k(phi) = k0 +k1*phi

#================= Explicit Euler =============================================
dphi = dt*(k0*np.dot(A,c) + k1*np.dot(np.diag(np.dot(B,c),k=0),np.dot(A,c)) + k1*np.dot(A1,c)**2)
dphi[0] = 0; dphi[-1] = 0
Bstar = np.copy(B)
Bstar[0] = 0; Bstar[-1] = 0

phi = np.zeros((len(t),len(xfine)))
for i in range(len(t)):
    interphi = np.dot(Bstar,c) + dphi + bc 
    c        = np.dot(Binv,interphi)
    dphi = dt*(k0*np.dot(A,c) + k1*np.dot(np.diag(np.dot(B,c),k=0),np.dot(A,c)) + k1*np.dot(A1,c)**2)
    phi[i]   = np.dot(C,c)
np.savetxt('phiKfunction',phi)

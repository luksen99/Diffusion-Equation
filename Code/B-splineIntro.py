
"""
Created on Fri Jun 29 10:49:45 2018
Introduction to basic B-splines
@author: Lukas
"""

import numpy as np
import matplotlib.pyplot as plt
from bspline import Bspline
from bspline.splinelab import aptknt, augknt

# Setup the spline

p      = 3       #order of the spline (should be at least 3 if 2nd order derivation is under investigation)
nknots = 10      # number of knots to generate (here endpoints count only once)
knots  = np.linspace(0,1,nknots) #create a knot grid with only physical points
k      = augknt(knots, p)        #add endpoints (ghostpoints) repeats as appropriate for spline order p
# alternatively: k = aptknt(knots,p)


Bspl = Bspline(k,p)
Bspl.plot()
plt.show()
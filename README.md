# Diffusion-Equation
In this Project we are providing a slover for Parabolic Partial Differential Equations, i.e. the heat diffusion equation, of the form: 


 ![Diffusion equation](https://latex.codecogs.com/gif.latex?C%28%5Cvarphi%28x%2Ct%29%29%5Ccdot%5Crho%5Cfrac%7B%5Cpartial%5Cvarphi%28x%2Ct%29%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%5Cleft%28k%28%5Cvarphi%28x%2Ct%29%29%5Ccdot%5Cfrac%7B%5Cpartial%5Cvarphi%28x%2Ct%29%7D%7B%5Cpartial%20x%7D%5Cright%29%20&plus;S%28t%2Cx%29.)
 
 Our approach is to use a combination of **finite elements** (B-Splines) to approximate the derivation in space and **explicit Euler** to approximate the evolution in time.
 
 ### Example showing our output:
 Temperature evolution in time and space after a laser pulse hits the probe
 
  Temperature evolution of probe |  Gaussian laser pulse S(x,t) hitting probe
:-------------------------:|:-------------------------:
 <img src="https://media.giphy.com/media/7TudjuaMsW2HP2xQ9Y/giphy.gif" width="320" height="300" />  |  <img src="https://github.com/luksen99/Diffusion-Equation/blob/master/Images/LaserPulse.png" width="320" height="300" />
 
 Here we consider a gaussian pulsed laser source S(x,t) = exp(-x)*G(t) hitting a probe in the middle. The probe gets heated up, as the pulse kicks in and the heat diffuses along the material until equilibrium is reached.


 
 

## Folders:
* 1) Documentation 
      
      
* 2) Descriptive \.ipynb files. To give a descriptive overview  of what can be solved with the package

 
* 3) Code \.py files containing code to comment, modify and work on

     
## Installation

use the download link

`python -m pip install --index-url https://test.pypi.org/simple/ solpde`

then run 

`from solpde import solpde`


### How to contribute: 
Fork from the `Developer`- branch and pull request to merge back into the original `Developer`- branch. 
Working updates and improvements will then be merged into the `Master` branch, which will always contain the latest working version.



With: 
* [Valentino Scalera](https://github.com/VaSca92)
* [UDCM Group of SU](http://udcm.fysik.su.se/)

Dependencies:

[Numpy](http://www.numpy.org/)

[Matplotlib](https://matplotlib.org/)

[Scipy.linalg](https://www.scipy.org/)

[B-splines](https://github.com/johntfoster)


`



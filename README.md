# Diffusion-Equation
In this Project we are providing a slover for Parabolic Partial Differential Equations, i.e. the heat diffusion equation, of the form: 


 ![Diffusion equation](https://latex.codecogs.com/gif.latex?C%28%5Cvarphi%28x%2Ct%29%29%5Ccdot%5Crho%5Cfrac%7B%5Cpartial%5Cvarphi%28x%2Ct%29%7D%7B%5Cpartial%20t%7D%20%3D%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20x%7D%5Cleft%28k%28%5Cvarphi%28x%2Ct%29%29%5Ccdot%5Cfrac%7B%5Cpartial%5Cvarphi%28x%2Ct%29%7D%7B%5Cpartial%20x%7D%5Cright%29%20&plus;S%28t%2Cx%29.)
 
 Our approach is to use a combination of **finite elements** (B-Splines) to approximate the derivation in space and **explicit Euler** to approximate the evolution in time.
 
 
 ![ALT TEXT](https://giphy.com/gifs/7TudjuaMsW2HP2xQ9Y/giphy.gif)

## Folders:
* 1) Descriptive \.ipynb files. To give a descriptive overview  of the problem and show intermediate results via Jupyter sessions

      *Recommeded to look at, for people being rather new to the project*
 
* 2) \.py files containing code to comment, modify and work on

     *Recommeded to look at, for people working on- and contributing to the code*
     
* 3) A report as a guide through the theoretical background of the code  

### How to contribute: 
Fork from the `Developer`- branch and pull request to merge back into the original `Developer`- branch. 
Working updates and improvements will then be merged into the `Master` branch, which will always contain the latest working version.



With: 
* [Valentino Scalera](https://github.com/VaSca92)
* [UDCM Group of SU](http://udcm.fysik.su.se/)



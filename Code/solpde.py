"""
Created on Mon Jul  9 13:51:45 2018

@author: Alber Lukas, Scalera Valentino 
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spl
from bspline import Bspline
from bspline.splinelab import aptknt
import time
from matplotlib.animation import FuncAnimation as movie
from mpl_toolkits.mplot3d import Axes3D





class heatpde(object): 
    
    def __init__(self): 
        self.xgrid          = 100           #number of points in x grid
        self.length         = 1             #length of x space,starting from 0
        self.start_time     = 0             #starting time (can be negative)
        self.final_time     = 3             #time when simulation stops
        self.time_step      = 0.0005        #timesteps in explicit Euler
        self.Left_BC_Type   = 0             #Boundary conditions Default is Neumann
        self.Right_BC_Type  = 0             #1=> Neumann; 0=> Dirichlet
        self.Left_BC        = lambda t: 0   #Value at Boundaries as fct of time
        self.Right_BC       = lambda t: 0
        self.init           = 0             # initial temperature of probe
        self.conductivity   =  lambda phi: 1          #Default for conductivity k(phi)
        self.heatCapacity   =  lambda phi: 1          #Default for het capacity C(phi)
        self.source         = 1     #source term for heating (can be a function of x and t)
        self.rho            = 1     #Density
        self.phi            = False #to not always recalculate
        
        
    def getProperties(self): # to depict the properties of the object
        for i in (self.__dict__): 
            print(i,' : ',self.__dict__[i])
        
        
    
    def Temp(self):         # load in all the defined properties
        num_of_points = 15  #Number of points used in the spline
        order         = 5   #order of the spline
        plt_points    = self.xgrid
        length        = self.length
        start_time    = self.start_time
        final_time    = self.final_time
        time_step     = self.time_step
        Left_BC_Type  = self.Left_BC_Type
        Right_BC_Type = self.Right_BC_Type
        Left_BC       = self.Left_BC
        Right_BC      = self.Right_BC
        conductivity  = self.conductivity 
        rfct          =  self.heatCapacity
        source        = self.source
        init          = self.init
        rho           = self.rho
        
        

        
        # Simulation Setup ====================================================
        # Capacity Function ---------------------------------------------------
        #dk/(dphi) is needed later in the core in every loop
        # evaluating the derivative of conductivity with respect to phi 
        # ---------------------------------------------------------------------
        
        def diff_conductivity(phi): 
            eps =1e-9
            dc = (conductivity(phi+eps)-conductivity(phi))/eps
            return(dc)
        
        # Capacity Function ---------------------------------------------------
        # evaluating coefficients of the passed on function via a system of linear eq.
        # ---------------------------------------------------------------------
        def capacity(r):
            A = np.array([[1,1,1,1],[1,2,4,8],[1,3,9,27],[1,4,16,64]])
            B = np.array([r(1),r(2),r(3),r(4)])
            rcoeff = np.linalg.solve(A,B)
            return(rcoeff)

        
        # Define Time and Space grid ------------------------------------------
        x = np.linspace(0, length , num_of_points)             # Space Grid
        t = np.arange(start_time, final_time, time_step)       #  Time Grid
        
        # Define Splines and Differentiation Matrices -------------------------
        knot_vector = aptknt(x, order)         # Creates knot points with Ghost points
        basis = Bspline(knot_vector, order)    # Generate a vector of Spline Objects
        A0 = basis.collmat(x, deriv_order = 0) # Generate Matrix A0 0st order derivative in space
        
        AA0 = basis.collmat(x, deriv_order = 0)# Generate a Boundary condition free A0 version
        AA0[-1,-1] = 1
        
        A1 = basis.collmat(x,deriv_order=1)    # Generate Matrix A1 1st order derivative in space
        A2 = basis.collmat(x, deriv_order = 2) # Generate Matrix A2 2st order derivative in space
        
        # Prepare "Smooth Plot Matrix"
        xx = np.linspace(0, length, plt_points)  # Gird to Plot
        C  = basis.collmat(xx)                   # Smooth Plot Matrix (read LaTeX notes)
        
        # Correct last spline  
        A0[-1,-1] = 1; C[-1,-1] = 1
        A1[-1]    = -np.flip(A1[0],0) # put first values of A1 to the last row
        
        # Prepare the inverse P to save time during the simulation ------------
        if  Left_BC_Type == 1 : A0[ 0] = A1[ 0] #set 1st and last row to 0
        if Right_BC_Type == 1 : A0[-1] = A1[-1] #needed to implement BC
        P = spl.inv(A0)
        
        # Modify Formulation to implement Boundary Condition ------------------
        A0[ 0] = 0;A1[ 0] = 0;A2[ 0] = 0 #First row is reserved for Boundary Condition
        A0[-1] = 0;A1[-1] = 0;A2[-1] = 0 # Last row is reserved for Boundary Condition
        
        # Time Evolution Matrix -----------------------------------------------
        M = np.dot( P, A0) + (time_step*np.dot( P, A2)) #only needed for the simple case
                                                        #see core
                                                        
        #initial c = coefficients for splines
        if isinstance(init,(int, float)) or len(init(x)) <len(x): #make the initial condition a function
            dummy = init                                          #in case input is a number 
            init  = lambda x: dummy + 0*x
        c = np.dot(spl.inv(AA0),init(x))

        # Prepare Boundary Condition according to which function is given to the class
        BC     = np.zeros((len(x),len(t)));
        BC[0]  = Left_BC(t)
        BC[-1] = Right_BC(t)
        #Prepare a matrix with the source data (space, time) to not always call 
        #the function in the loop, see core
        if isinstance(source,(int,float)): # make the source a function
            dummy1 = source                #in case the input is a number, i.e.a constant
            source = lambda x,t: dummy1 + 0*x +0*t
            
        xmg,tmg = np.meshgrid(x,t)
        sourceM = source(xmg,tmg)
        sourceM[:,0] = 0; sourceM[:,-1] = 0 #set last and first row 0 for BC
        
        #Prepare Array to store results
        phi = np.zeros((len(t),len(xx)))
        
        # End of Simulation Setup =============================================
        
        # MAIN LOOP -----------------------------------------------------------
        # =====================================================================
        # Decide which case is relevant and solve the according for loops in core.py file
        # Depending on which _BC_Type (either Neumann or Dirichlet) is radsed, 
        #respective boundary conditions are taken into consideration. 
        # =====================================================================

        if Left_BC_Type == 0 and Right_BC_Type == 0: #Dirichlet on both sides
            print('Dirichlet condition on both sides')
            if isinstance(conductivity, (int, float)) == True and isinstance(rfct, (int, float)) == True: 
            #k(phi) = k0 & r(phi) = r0 Conditions check if conductivity & capacity are constants
                print('Constant Conductivity and Capacity and Dirichlet boundary conditions')
                print('No source and no density is taken under consideration')
                k0  = conductivity; r0 = rfct
                phi = core.simple_DD(M,t,c,k0,P,BC,C,phi)
            if isinstance(conductivity, (int, float)) == False and isinstance(rfct, (int, float)) == True:    
            #Conductivity: k(phi) = k0 + k1*phi and Capacity: r(phi) = r0
                print(r'Generic $k(\phi)$ and Capacity:$r(\phi) = r0$')
                r0  = rfct
                phi = core.genK_phi_DD(t,c,A0,A1,A2,diff_conductivity,conductivity,sourceM,r0,time_step,P,C,BC,phi,rho,x)
        
            if isinstance(conductivity, (int, float)) == False and isinstance(rfct, (int, float)) == False:
            #Conductivity: k(phi) & Capacity: r(phi) are both generic
                print(r'Generic Conductivity: $k(\phi)$ and Capacity: $r(\phi)$')
                r0  = capacity(rfct)[0]; r1 = capacity(rfct)[1]; r2 = capacity(rfct)[2]; r3 = capacity(rfct)[3]
                phi = core.genK_genR_DD(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,rho,x)
        
               
                
        if Left_BC_Type == 1 and Right_BC_Type == 0:# Neumann condition on RHS: 
            print('Left side: Neumann- ; Right side: Dirichlet BC')
            print(r'Generic Conductivity: $k(\phi)$ and Capacity: $r(\phi)$')
            side = 0
            r0   = capacity(rfct)[0]; r1 = capacity(rfct)[1];r2 = capacity(rfct)[2]; r3 = capacity(rfct)[3]
            phi  = core.genK_genR_ND(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,side,rho,x)
            
        if Left_BC_Type == 0 and Right_BC_Type == 1:# Neumann condition on RHS: 
            print('Left side: Dirichlet- ; Right side: Neumann Boundary condition')
            print(r'Generic Conductivity: $k(\phi)$ and Capacity: $r(\phi)$')
            side = -1
            r0   = capacity(rfct)[0]; r1 = capacity(rfct)[1];r2 = capacity(rfct)[2]; r3 = capacity(rfct)[3] 
            phi  = core.genK_genR_DN(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,side,rho,x)                    
                
        if Left_BC_Type == 1 and Right_BC_Type == 1:# Neumann condition on RHS & LHS: 
            print('Both sides Neumann Boundary condition')
            print(r'Generic Conductivity: $k(\phi)$ and Capacity: $r(\phi)$')
            r0  = capacity(rfct)[0]; r1 = capacity(rfct)[1];r2 = capacity(rfct)[2]; r3 = capacity(rfct)[3]
            phi = core.genK_genR_NN(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,rho,x) 
         
        return(phi) # in every case phi gets returned by the core
                    # the values of phi are stored in self.Temp as a Matrix
                    
                    
    # Visualization -----------------------------------------------------------
    #in order to visualize there are a few ready made visualization functions
    #therefor the class 'visual' gets called and the properties of consideration
    #defined before are arguments of the functions
    #--------------------------------------------------------------------------
    
    def spaceVStime(self): #2D countour plot
        xx = np.linspace(0,self.length,self.xgrid)
        t  = np.arange(self.start_time, self.final_time, self.time_step)
        if isinstance(self.phi,bool): #to not calculate again, once propertie is obtained
            self.phi = self.Temp()
        visual.spaceVStime(xx,t,self.phi)

    def single_point(self,myspace): #one point which can be selectet vs entire time
        t  = np.arange(self.start_time, self.final_time, self.time_step)
        xx = np.linspace(0,self.length,self.xgrid)
        k  = np.where(np.min(xx-myspace))[0][0]
        if isinstance(self.phi,bool): #to not calculate again, once propertie is obtained
            self.phi = self.Temp()
        visual.single_point(t,self.phi,k,myspace) 
        
    def single_time(self,mytime): #snapshot of one specific time vs entire x-space
        xx = np.linspace(0,self.length,self.xgrid)
        t  = np.arange(self.start_time, self.final_time, self.time_step)
        k  = np.where(np.min(t-mytime))[0][0]
        if isinstance(self.phi,bool): #to not calculate again, once propertie is obtained
            self.phi = self.Temp()
        visual.single_time(xx,self.phi,k,mytime)         

    def animated(self,speed): #animation showing the evolution in x-space as time passe by
        xx = np.linspace(0,self.length,self.xgrid)
        t  = np.arange(self.start_time, self.final_time, self.time_step)
        if len(t) % speed != 0:
            speed = speed-(len(t)%speed)
        if isinstance(self.phi,bool): 
            self.phi = self.Temp()        
        visual.animated(xx,t,self.phi,speed)
        
       
    def surface(self): #surface plot, similar to spaceVStime just in 3D
        xx = np.linspace(0,self.length,self.xgrid)
        t  = np.arange(self.start_time, self.final_time, self.time_step)
        if isinstance(self.phi,bool): 
            self.phi = self.Temp()        
        visual.surface(xx,t,self.phi)  
        
    def surf_source(self): #3D surface plot of the input source
       x = np.linspace(0,self.length,50)
       t = np.linspace(self.start_time, self.final_time, 50)               
       xmg, tmg = np.meshgrid(x, t)
       source = self.source
       sourceM = source(xmg,tmg)
       visual.vis_source(xmg,tmg,sourceM)     
        
   
        
        
        # END Visualization ----------------------------------------------------
        
        
# =============================================================================
# Solvers for the loops accoridng to which boundary condition is set        
# =============================================================================
        
class core(object):
# =====================================================================
# Simulation CORE
# here all the loops and the entire evolution in time is calculated
# description of the loop is in genK_phi_DD()- function    
# =====================================================================
            
    def simple_DD(M,t,c,k0,P,BC,C,phi): 
        """
        Most simple case, where the BC_Type of both sides is Dirichlet
        and there is no Source, no density, no capacity C(phi)
        """
        start = time.time()
        for i in range(len(t)):
            c      = k0*np.dot( M, c) + np.dot( P, BC[:,i])
            phi[i] = np.dot(C,c)
        end = time.time()
        print('Elapsed time: ', end-start)     
        return(phi) 
    
    def  genK_phi_DD(t,c,A0,A1,A2,diff_conductivity,conductivity,sourceM,r0,time_step,P,C,BC,phi,rho,x):
        """
        Here a generic k(phi), i.e. conductivity is considered. heatCapacyty = constant = r0
        Boundary conditions are of Dirichlet type on both sides
        """
        start = time.time()
        for i in range(len(t)): 
            phi0 = np.dot(A0,c); phi1 = np.dot(A1,c); phi2 = np.dot(A2,c) #phi as dot product of spline matrix and coefficents
            Flow_1 = diff_conductivity(phi0) * phi1**2                    #flow_1 (see theory paper)
            Flow_2 = conductivity(phi0) * phi2                            #flow_2 (see theory paper)
            dphi   = (Flow_1 + Flow_2 + sourceM[i])/(r0*rho)              #partial derivative of phi
            intphi = phi0 + time_step*dphi + BC[:,i]                      #intermediate phi 
            c      = np.dot(P,intphi)                                     #multiply with the inverse P from the left side to obtain coefficients c
            phi[i] = np.dot(C,c)                                          #map on a more fine Spline Matrix C, for more smooth plotting
        end = time.time()
        print('Elapsed time: ', end-start)       
        return(phi)
    
    def genK_genR_DD(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,rho,x):
        """
        General case of Dirichlet boundary conditions on both sides
        """
        start = time.time()
        for i in range(len(t)):
            phi0   = np.dot(A0,c); phi1 = np.dot(A1,c); phi2 = np.dot(A2,c)
            Flow_1 = diff_conductivity(phi0) * phi1**2 
            Flow_2 = conductivity(phi0) * phi2
            r_phi  = 1/(r0 + r1*np.dot(AA0,c) + r2*np.dot(AA0,c)**2 + r3*np.dot(AA0,c)**3)
            dphi   = r_phi*(Flow_1 + Flow_2 +  sourceM[i])/rho
            intphi = phi0 + time_step*dphi + BC[:,i] 
            c      = np.dot(P,intphi)
            phi[i] = np.dot(C,c)
        end = time.time()
        print('Elapsed time: ', end-start) 
        return(phi)   
        
    def genK_genR_ND(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,side,rho,x):
        """
        General case ofNeumann BC LHS Dirichlet boundary conditions on RHS
        Therefor the Boundary conditions on side = 0 => lhs need to get 
        updated in the loop.
        """
        start = time.time()
        for i in range(len(t)):
            phi0 = np.dot(A0,c); phi1 = np.dot(A1,c); phi2 = np.dot(A2,c)
            Flow_1  = diff_conductivity(phi0) * phi1**2
            Flow_2  = conductivity(phi0) * phi2    
            r_phi   = 1/(r0 + r1*np.dot(AA0,c) + r2*np.dot(AA0,c)**2 + r3*np.dot(AA0,c)**3)
            dphi    = r_phi*(Flow_1 + Flow_2 +  sourceM[i])/rho 
            BC[side,i] = BC[side,i] /conductivity(c[side]) #updating the boundary condition 
            intphi  = phi0 + time_step*dphi + BC[:,i]      #to the Neumann case 
            c       = np.dot(P,intphi)
            phi[i]  = np.dot(C,c)    
        end = time.time()
        print('Elapsed time: ', end-start) 
        return(phi) 
        
    def genK_genR_DN(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,side,rho,x):
        """
        General case ofNeumann BC RHS Dirichlet boundary conditions on LHS
        Therefor the Boundary conditions on side = -1 => rhs need to get 
        updated in the loop.
        """
        start = time.time()
        for i in range(len(t)):
            phi0 = np.dot(A0,c); phi1 = np.dot(A1,c); phi2 = np.dot(A2,c)
            Flow_1  = diff_conductivity(phi0) * phi1**2
            Flow_2  = conductivity(phi0) * phi2    
            r_phi   = 1/(r0 + r1*np.dot(AA0,c) + r2*np.dot(AA0,c)**2 + r3*np.dot(AA0,c)**3)
            dphi    = r_phi*(Flow_1 + Flow_2 +  sourceM[i])/rho 
            BC[side,i] = BC[side,i] /conductivity(c[side]) 
            intphi  = phi0 + time_step*dphi + BC[:,i] 
            c       = np.dot(P,intphi)
            phi[i]  = np.dot(C,c)    
        end = time.time()
        print('Elapsed time: ', end-start) 
        return(phi) 
    
    
    def genK_genR_NN(t,c,A0,AA0,A1,A2,diff_conductivity,conductivity,sourceM,r0,r1,r2,r3,time_step,P,C,BC,phi,rho,x):
        """
        General case ofNeumann BC on both sides
        Therefor the Boundary conditions need to get updated in the loop.
        """        
        lhs = 0; rhs = -1
        start = time.time()
        for i in range(len(t)):
            phi0 = np.dot(A0,c); phi1 = np.dot(A1,c); phi2 = np.dot(A2,c)
            Flow_1  = diff_conductivity(phi0) * phi1**2
            Flow_2  = conductivity(phi0) * phi2    
            r_phi   = 1/(r0 + r1*np.dot(AA0,c) + r2*np.dot(AA0,c)**2 + r3*np.dot(AA0,c)**3)
            dphi    = r_phi*(Flow_1 + Flow_2 +  sourceM[i])/rho 
            BC[lhs,i] = BC[lhs,i] /conductivity(c[lhs]) 
            BC[rhs,i] = BC[rhs,i] /conductivity(c[rhs]) 
            intphi  = phi0 + time_step*dphi + BC[:,i] 
            c       = np.dot(P,intphi)
            phi[i]  = np.dot(C,c)  
        end = time.time()
        print('Elapsed time: ', end-start) 
        return(phi) 
        
# =============================================================================
# End of solvers        
# =============================================================================


#Visuals ======================================================================
class visual(object):   
    # =========================================================================
    # 2D countour plot where x = space, y= time and z = temperature
    # =========================================================================
    def spaceVStime(x,t,phi):
        plt.contourf(x,t,phi,50,cmap = 'plasma')
        plt.colorbar( orientation='vertical', shrink=0.8)
        plt.xlabel('x-Space')
        plt.ylabel('time')
        plt.title(r'$\phi(t,x)$')
        plt.show()
    
    # =========================================================================
    # One point which can be selectet vs entire time. 
    # x = time; y = temperature at selected point in x space
    # =========================================================================
    def single_point(t,phi,i,mypoint):
        plt.plot(t,phi[:,i])
        plt.ylabel('Amplitude')
        plt.xlabel('Time')
        plt.title('Evolution of Amplitude at one point x='+'{:10.4}'.format(mypoint)+' m in space')
        plt.grid()
        plt.show()
        
    # =========================================================================
    # Snapshot of one specific time vs entire x-space
    # x = space, y = temperature at selected time i
    # =========================================================================
    def single_time(x,phi,i,mytime):
        plt.plot(x,phi[i,:])
        plt.ylabel('Amplitude')
        plt.xlabel('Space')
        plt.title('Evolution of Amplitude at time ='+'{:10.4}'.format(mytime)+'s')
        plt.grid()
        plt.show()        
    
    # =========================================================================
    # Animation showing the evolution in x-space as time passe by
    # here we use: from matplotlib.animation import FuncAnimation as movie
    # =========================================================================
    def animated(x,t,phi,speed):
        fig, ax = plt.subplots()
        line,   = ax.plot([], [], 'r', animated=True)
        ax.set_xlim(0, x[-1]); ax.set_ylim(np.min(phi)-(np.mean(phi)-np.min(phi))/2, np.max(phi)+(np.max(phi)-np.mean(phi))/2)
        time_text = ax.text(0.02, 0.95, "", transform = ax.transAxes)
        plt.xlabel('Depth of Material'); plt.ylabel('Temperature')
        plt.title('Evolution of Temperature in space and time')
    
        def update(frame):
            line.set_data(x,phi[speed*frame])
            time_text.set_text("time = " + str(t[speed*frame]) + " s")
            return line, time_text
    
        ani = movie(fig, update, blit=True)
        plt.grid()
        plt.show()  
 
    # =========================================================================
    # surface plot, similar to spaceVStime just in 3D
    # here we use: from mpl_toolkits.mplot3d import Axes3D
    # =========================================================================       
    def surface(x,t,phi):
        x, t = np.meshgrid(x[::2], t[::10])         
        fig  = plt.figure()
        ax   = fig.gca(projection='3d')
        surf = ax.plot_surface(x,t,phi[::10,::2],cmap = 'plasma')
        fig.colorbar(surf,shrink=0.7, aspect=5)
        plt.xlabel('x-Space')
        plt.ylabel('time')
        plt.title(r'$\phi(t,x)$')
        plt.show()
        
    # =========================================================================
    # 3D surface plot of the input source
    # here we use: from mpl_toolkits.mplot3d import Axes3D
    # =========================================================================         
    def vis_source(x,t,s): 
        fig  = plt.figure()
        ax   = fig.gca(projection='3d')
        surf = ax.plot_surface(x,t,s,cmap = 'jet')
        fig.colorbar(surf,shrink=0.7, aspect=5)
        plt.xlabel('x-Space')
        plt.ylabel('time')
        plt.title(r'S(x,t)')
        plt.show()    



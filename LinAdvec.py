import numpy as np
import matplotlib.pyplot as plt


#function defining the intial and analytic solution
def initialBell(x):
    return np.where(x%1. < 0.5, np.power(np.sin(2*x*np.pi),2),0)
    #%n remainder after division by n

#put everything inside a main function to avoid global variables    
def main():
    # Setup space, initial phi profile and Courant num
    nx = 40             # number of points in space
   # c = 0.2             #Courant number
        #spatial variable going from zero to one inclusive
    x = np.linspace(0.0,1.0,nx+1)
    
    arraysize = 100
    #dt = 0.01
    carray = np.zeros([arraysize])
    delta = np.zeros([arraysize])
    dtlist = np.linspace(0.001,0.01,arraysize)
    for h in range(0,arraysize):
        
        dt = dtlist[h]
        # Set up changes
        
        #given time
        t = 10
        
        
       # nt = 40
        nt = round(t/dt)
        #derived quantities
        u = 1.
        dx = 1./nx
        c = u * dt/dx
       # dx = u*dt./c
        
        #dt = c*dx/u
       # t = nt*dt
        
        
        

        #Three time levels of the dependent variable, phi
        phi = initialBell(x)
        phiNew = phi.copy()
        phiOld = phi.copy()
        
        #FTCS for the first time-step, looping over space
        
        for j in range(1,nx):
            phi[j] = phiOld[j] - 0.5*c*(phiOld[j+1] - phiOld[j-1])
        #Apply periodic boundary conditions
        phi[0] = phiOld[0] - 0.5*c*(phiOld[1] - phiOld[nx - 1])
        phi[nx] = phi[0]
        
        #loop over the remaining time-steps (nt) using CTCS
        
        
        
        for n in range(1,int(nt)):
        #loop over space
            for j in range(1,int(nx)):
                phiNew[j]= phiOld[j] - c*(phi[j+1] - phi[j-1])
             #apply periodic boundary conditions
            phiNew[0] = phiOld[0] - c*(phi[1] - phi[nx-1])
            phiNew[nx] = phiNew[0]
             #update phi for the next time-step
            phiOld = phi.copy()
            phi = phiNew.copy()
             
      
         
        #plot the solution in comparison to the analytic solution
        plt.plot(x, initialBell(x - u*t), 'k', label = 'analytic')
        plt.plot(x,phi,'b',label = 'CTCS')
        plt.legend(loc = 'best')
        plt.ylabel('$\phi$')
        plt.axhline(0,linestyle = ':',color = 'black')

     
     
     #calculate the difference between the plots
        diff = abs(sum(initialBell(x - u*t) - phi)) 
        delta[h] = diff
        carray[h] = c
        print('Difference between is ', diff,' with c = ',c)
 
    #plt.show()
main()
         
         
         
         
         

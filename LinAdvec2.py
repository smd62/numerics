import numpy as np
import matplotlib.pyplot as plt


#function defining the intial and analytic solution
def initialBell(x):
    return np.where(x%1. < 0.5, np.power(np.sin(2*x*np.pi),2),0)
    #%n remainder after division by n

#put everything inside a main function to avoid global variables    
def main():
    # Setup space, initial phi profile and Courant num
    arraysize = 100
    nxArray = np.linspace(10,200,arraysize)
                   #spatial variable going from zero to one inclusive
    
    
    ntArray = np.linspace(10,1000,arraysize)
      #given time
    t = 10.
    #nx = 40             # number of points in space
   # c = 0.2             #Courant number


    
   
    #dt = 0.01
    carray = np.zeros([arraysize,arraysize])
    delta = np.zeros([arraysize,arraysize])
    dtlist = np.linspace(0.001,0.025,arraysize)
    for tindex in range(0,arraysize):
        for xindex in range(0,arraysize):
            
            nx = int(nxArray[xindex])
            nt = int(ntArray[tindex])
            
            dt = t/nt
            x = np.linspace(0.0,1.0,nx+1)
            # Set up changes
            
          
            
            
           # nt = 40
            #nt = round(t/dt)
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
                phiNew[0] = phiOld[0] - c*(phi[1] - phi[int(nx)-1])
                phiNew[int(nx)] = phiNew[0]
                 #update phi for the next time-step
                phiOld = phi.copy()
                phi = phiNew.copy()
                 
          
             
            #plot the solution in comparison to the analytic solution
           # plt.plot(x, initialBell(x - u*t), 'k', label = 'analytic')
            #plt.plot(x,phi,'b',label = 'CTCS')
            #plt.legend(loc = 'best')
            #plt.ylabel('$\phi$')
            #plt.axhline(0,linestyle = ':',color = 'black')

         
         
         #calculate the difference between the plots
            diff = abs(sum(initialBell(x - u*t) - phi)) 
            
            
            if c > 1:
               diff = np.nan
               
                
            
            delta[xindex,tindex] = diff
            carray[xindex,tindex] = c
          #  print('dx = ',dx,' with dt = ',dt)
           # print('Difference between is ', diff,' with c = ',c)
    print('delta = ',delta)
    #print('c = ', carray)
    #print('nx = ',nxArray)
    #print('nt = ',ntArray)
    
    print('len(nx) = ', len(nxArray))
    print('len(nt) = ', len(ntArray))
    print('len(delta) = ', np.shape(delta))


    X, Y = np.meshgrid(nxArray, ntArray)
    
    plt.figure()
    cp = plt.contourf(X, Y, delta, 100,cmap = 'Reds',vmax = 2)
    cbar = plt.colorbar(cp)
    cbar.ax.set_ylabel('Difference')
    plt.xlabel('nx')
    plt.ylabel('nt')
    plt.title('Difference between analytical and CTCS')
    plt.show()
    
    
    
    
    plt.figure()
    cp2 = plt.contourf(X, Y, carray, 100,cmap = 'gist_earth_r')
    cbar = plt.colorbar(cp2)
    cbar.ax.set_ylabel('c')
    plt.xlabel('nx')
    plt.ylabel('nt')
    plt.title('Courant Number')
    plt.show()    
    
    
   # plt.pcolor([nxArray, ntArray],delta)
   # plt.xlabel(
   # tick_spacing = 10
   # plt.yticks(np.arange(0, len(ntArray), tick_spacing), ntArray[0::tick_spacing])
   # plt.xticks(np.arange(0, len(nxArray), tick_spacing), nxArray[0::tick_spacing])   
   
  #  plt.show()
    
    
    #plt.show()
main()
         
         
         
         
         

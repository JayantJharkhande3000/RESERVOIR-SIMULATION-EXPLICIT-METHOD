# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 13:08:47 2020

@author: Dell
"""

import numpy as np
import matplotlib.pyplot as plt 

print(   "\t ########################################################################### \t"  )
print(   "\t ############### Reservoir Simulation Explicit solution to 1D  ###################### \t"  )
print(   "\t ######## For dirichlet(constant pressure) boundary consition   ############ \t"  )
print(   "\t ############ left side constant pressure boundary consition ############### \t"  )
print(   "\t ################## Right side no flow  boundary consition   ############### \t"  )
print(   "\t ########################################################################### \t"  )
    



L = 10000
print("\n\tThe lenght of the reservoir is ", str(L))
n = 4 
print("\n\tBlock node in the reservoir is ", str(n))
P0 = 1000
print("\n\tThe intial pressure of the reservoir is ", str(P0))
P_left = 2000
print("\n\tThe pressure at the left boundary of the reservoir is ", str(P_left))
dx = L/n
print("\n\tThe lenght of Blocks in the reservoir is ", str(dx))

porosity = 0.2
print("\n\tthe porosity value of the reservoir is ", str(porosity))

permeability  = 50
print("\n\tthe permeability(md) value of the reservoir is ", str(permeability))

viscosity = 1
print("\n\tthe viscosity value is ", str(viscosity))

area = 200000
print("\n\tCross sectional area of the reservoir ", str(area))

compressibility = 1*10**(-6)
print("\n\tcompressibility of the reservoir is ", str(compressibility))


#################### final time for simulation is  ##############################
t_final = 3
print("\n\t the reservoir simulation time in days is  ", str(t_final))

#################### time step  ###############################################
dt = 1
print("\n\t the reservoir simulation incremental time step in days is "+ str(dt)+ "day")

########  x represent the lenght distribution array in the L direction with n nodes ########
x = np.linspace(dx/2, L-dx/2, n)


########  P represent the pressure distribution array in the reservoir with n nodes ########
P = np.ones(n)*P0

###### dPdt represent the pressure_gradient update distribution array in the reservoir with n nodes ########
dPdt = np.empty(n)

############  T represent the time distribution array in with the dt increment #############
t =  np.arange(0 , t_final, dt)



alpha = float(permeability /(porosity*viscosity*compressibility))

print("\n\t the neta constant value is  ", str(alpha))

neta = float((alpha*dt)/(dx**2))
print("\n\t the neta constant value is  ", str(neta))



neta = 0.2532
print("\n\t the neta constant value is  ", str(neta))


print("\n############## pressure distribution ################\n")
print("pressure distribution at day 0 is ", str(P))
    
for j in range(1, len(t)):
#    plt.clf()
    Pin = 2*P_left - P[0]
    dPdt[0] = (neta*(Pin - 2*P[0] + P[1] ))
    #print("\nblock = 0th pressure gradient update = ",dPdt[0])
    #print("Pin value = ",Pin)
    for i in range(1,n-1):
                        
        dPdt[i] = (neta*(P[i-1] - 2*P[i] + P[i+1] ))

    dPdt[n-1] = (neta*(P[n-2] - 2*P[n-1] + P[n-1] ))
    #print("block = n-1th pressure gradient update = ",dPdt[n-1])
    P = P + dPdt
    if j <= 10:
        print("pressure distribution at day "+str(j)+ " is =",P)
    elif j>=280:
        print("pressure distribution at day "+str(j)+ " is =",P)
plt.figure(1)
plt.plot(x,P)
plt.show
#    

#    plt.figure(1)
#    plt.plot(x,P)
#    plt.show
#    plt.pause(0.1)
#        
   
#print(P) 

    
    
    
    
    
    
    
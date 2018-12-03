"""
Created on Sun Nov  4 16:36:45 2018

@author: jeremyuzan
"""

import numpy as np
import matplotlib.pyplot as plt



List_NX = [4, 5, 7, 8, 10, 15, 25]

NX = 20
CFL = 0.450
dx = 1/NX
dt = CFL * dx**2
T = 0.4
alpha = 1
beta = 1
U_data = [alpha * np.sin(2*np.pi*i*dx) + beta * np.cos(2*np.pi*i*dx) for i in range(NX+1)]
U_sol = [np.exp(-4*np.pi*np.pi*T)*U for U in U_data]





def scheme_1(dt,dx,CFL,U_old,NX):# gives dU
    ddU = np.zeros(NX+1)
    if(np.abs(CFL-(dt/(dx*dx)))>10**(-10)):
        print("la condition CFL ne fonctionne pas")
        return ddU
    
    for j in np.arange(2, NX-1):

        ddU[j] =(4/3)*CFL*(U_old[j+1] - 2*U_old[j]+ U_old[j-1])-(1/12)*CFL*(U_old[j+2] - 2*U_old[j] + U_old[j-2])+0.5*CFL*CFL*(U_old[j+2] - 4*U_old[j+1] + 6*U_old[j] - 4*U_old[j-1] + U_old[j-2])
        
       
    ddU[NX-1] = CFL*4/3*(U_old[0] - 2*U_old[NX-1]+ U_old[NX-2])+CFL*(-1/12)*(U_old[1] - 2*U_old[NX-1] + U_old[NX-3])+CFL*CFL*0.5*(U_old[1] - 4*U_old[0] + 6*U_old[NX-1] - 4*U_old[NX-2] + U_old[NX-3])

    ddU[1] = CFL*4/3 * (U_old[2] - 2 * U_old[1] + U_old[0])+CFL*(- 1/12) * (U_old[3] - 2 * U_old[1] + U_old[NX-1])+CFL*CFL*0.5*(U_old[3] - 4*U_old[2] + 6*U_old[1] - 4*U_old[0] + U_old[NX-1])
    
    ddU[0] = CFL*4/3 * (U_old[1] - 2 * U_old[0] + U_old[NX-1])+CFL*(- 1/12) * (U_old[2] - 2 * U_old[0] + U_old[NX-2])+CFL*CFL*0.5*(U_old[2] - 4 * U_old[1] + 6 * U_old[0] - 4 * U_old[NX-1] + U_old[NX-2])

    return ddU






def rec_scheme(dt,dx,U_data,NX,T,scheme):
    t = 0
    U = U_data.copy()
    CFL = (dt/(dx*dx))
    n = 0
    N = int(T/(10*dt))
    
    while(t<T):
        U += scheme(dt,dx,CFL,U,NX)
        U[NX] = U[0] 
        t += dt
        n+=1
        if(n%N == 0):
            print("t = ",t)
    
    return U





def L2_diff(dx,NX,U_1,U_2):
    L22 = 0
    
    for i in range(NX):
        L22 += dx*((U_1[i]-U_2[i])**2)
    
    return np.sqrt(L22)






Err_list = []
U_data = [alpha * np.sin(2*np.pi*i*dx) + beta * np.cos(2*np.pi*i*dx) for i in range(321)]
U_sol = [np.exp(-4*np.pi*np.pi*T)*U for U in U_data]
L2_std = L2_diff(1/320,320,U_sol,np.zeros(321))

for NX in List_NX: 
    dx = 1/NX
    dt = CFL * dx**2
    
    print("NOMBRE D'OPERATIONs : ", int(T/dt)," schemes computations")
    U_data = [alpha * np.sin(2*np.pi*i*dx) + beta * np.cos(2*np.pi*i*dx) for i in range(NX+1)]
    U_sol = [np.exp(-4*np.pi*np.pi*T)*U for U in U_data]
    U_num = rec_scheme(dt,dx,U_data,NX,T,scheme_1)
    err = L2_diff(dx,NX,U_sol,U_num)/L2_std
    Err_list.append(err)

    print("POUR ",NX,"POINTS DE GRILLE, L'ERREUR DE NORME L^2 EST ",err)
    plt.plot(U_num,"r", label = "SOLUTION NUMERIQUE")
    plt.plot(U_sol,"b", label = "SOLUTION THEORIQUE")
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    
    plt.show()




#Calcul de l'erreur numérique en fonction du nombre de pas de grille
print("Erreur en norme L^2:                 ",Err_list)
inversef = [1/x for x in List_NX]
#plt.figure(0)


plt.plot(List_NX, inversef, "black", label = "Erreur L^2 f(x)=1/x")
plt.plot(List_NX, Err_list, "yellow", label = "Err L^2 théo VS num ")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0.)









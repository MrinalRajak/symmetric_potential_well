# Symmetric potential weel eigenvalues and eigenfunctions

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import simps
from scipy.optimize import brentq
a=4*(1./0.529177) # width of the potential well(here well width is 4 angstrom)
N=10000 # number of points to take
Vo=(14./27.211) # parameter determining depth of the well(14 electron volt depth assumed)
en=np.linspace(0.00000001,Vo,100000) # vector of energies where we look for the stable states
m=1
h_bar=1
v=np.sqrt(2*m*en/(h_bar**2))*(a/2.)
uo=np.sqrt(m*a**2*Vo/(2.*h_bar**2))
f_sym=lambda v:v*np.tan(v)-np.sqrt((uo)**2-(v)**2) # symmetric case
f_asym=lambda v: (-1./np.tan(v))*v-np.sqrt((uo)**2-(v)**2) # Antisymmetric case
def f_sym1(v):
    return v*np.tan(v)   # symmetric case1
def f_circle(v):
    return np.sqrt((uo)**2-(v)**2) # the equation of circle
def f_asym1(v):
    return (-1./np.tan(v))*v   # antisymmetric case1
plt.grid(True)
plt.plot(v,f_sym1(v),v,f_asym1(v),v,f_circle(v))
plt.ylim(0,12.1)
plt.show()

v_zeros_s=[]   # stores zeros of v symmetric
v_zeros_a=[]   # stores zeros of v antisymmetric
v_zeros=[]     # stores zeros of v
Eigenenergies=[]

# Find the zeros for the symmetrical case
s=np.sign(f_sym(v))
for i in range(len(s)-1):
    if s[i]+s[i+1]==0:
        zero=brentq(f_sym,v[i],v[i+1])
        v_zeros_s.append(zero)
        l1=v_zeros_s[::2]  #selecting non divergent points of tan(v)
# Find the zeros for the antisymmetrical case
s=np.sign(f_asym(v))
for i in range(len(s)-1):
    if s[i]+s[i+1]==0:
        zero=brentq(f_asym,v[i],v[i+1])
        v_zeros_a.append(zero)
        l2=v_zeros_a[::2] # selecting non divergent points of tan(v)
v_zeros=l1+l2 # found zeros
print(v_zeros)

for i in v_zeros:
    Eigenenergies.append(2*i**2*h_bar**2/(m*a**2))
print("Atomic units",Eigenenergies)
Eigenenergies=np.asarray(Eigenenergies)
print("Electron Volts",27.211*Eigenenergies)

#E=Eigenenergies[0] # selecting eigenenergy for the ground state(Even/symmetric)
E=Eigenenergies[1]  # selecting eigenenergy for the 1st excited state(Even/symmetric)
E=Eigenenergies[2]  # selecting eigenenergy for the 1st excited state(Even/antisymmetric)

'''
potential function in the finite square well.width is a and depth is global variable Vo
'''
def V(x):
    a=4*(1./0.529177)
    if abs(x)>a/2.:
        return Vo
    else:
        return 0.
def schrodinger(w,x):
    psi=w[0]
    psi1=w[1]
    dpsidx=psi1
    dpsi1dx=-2*m*((E-V(x))/(h_bar**2))*psi
    return np.array([dpsidx,dpsi1dx])
#w=[1.,0.] # Initial condition for the Even states/symmetric states
w=[0.,1.] # Initial condition for the odd states/asymmetric states
xp=np.linspace(0,a,N) # positive x-axis
sol=odeint(schrodinger,w,xp) #(function,initial y,input x array)
psi_plus=sol[:,0]
psi1_plus=sol[:,1]
xn=np.linspace(0,-a,N) # negative x-axis
sol1=odeint(schrodinger,w,xn) # (function,initial y,input x array)
psi_minus=sol[:,0]
psi1_minus=sol[:,1]
# Construction of the final psix and psi
psix=np.concatenate((xn[::-1],xp[1:]))
psi=np.concatenate((psi_minus[::-1],psi_plus[1:]))
# normalization constant using simpson rule from scipy
sq_arr_psi=psi**2
print("scipy based simpson integration result as normalization constant=",simps(sq_arr_psi,dx=(xp[1]-xp[0])))
N=1./simps(sq_arr_psi,dx=(xp[1]-xp[0]))
# normalization wavefunction
psi=N*psi
plt.plot(psix,psi,'o',c='0.3')
plt.show()
                   




































































































        
        







































        































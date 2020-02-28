#Module for the 2D potential. Potential modes are:
#1. 2D Harmonic oscillator 
#2. 2D finite well
#3. 2D Gaussian 
#4. 2D Morse or anharmonic potential
import numpy as np
from numpy import pi,sqrt
mass=1.
###########################################################################
#1. 
k_x=1.
k_y=1.
def potential_ho(x,y):
    v=.5*k_x*x*x+.5*k_y*y*y+15*np.exp(-(x*x+y*y))
#    +15*np.exp(-(x*x+y*y))
    return v
###########################################################################
#2.
#sides of the well
ax=2.
ay=2.
#ground state energy
E0=pi*pi*(1/ax**2+1/ay**2)/2
#height of the well
v0=5.*E0
def potential_fwell(x,y):
    if x<=-ax/2 and y<=-ay/2:
        v=v0
    elif x<=-ax/2 and y>=ay/2:
        v=v0
    elif x>=ax/2 and y<=-ay/2:
        v=v0
    elif x>=ax/2 and y>=ay/2:
        v=v0
    else:
        v=0.
    return v
    
        
    
    
#define a function to compute the normalized Harmonic oscillator wavefunctions
#size of the Harmonic function basis in each dimension
import numpy as np
from numpy import pi,sqrt
import pot2D
sq2=sqrt(2.)
#x parameters
omega_x=1.
m_x=1.
alpha_x=m_x*omega_x
#y parameters
omega_y=1.
m_y=1.
alpha_y=m_y*omega_y
#Estimate maximum lengths along each dimension for a given basis size
def xy_dim(nsize_x,nsize_y):
    xmax=sqrt(2.*(nsize_x+.5)/(m_x*omega_x))*1.5
    ymax=sqrt(2.*(nsize_y+.5)/(m_y*omega_y))*1.5
    return xmax,ymax

def harmonicfun_x(x,n):
    chi=sqrt(alpha_x)*x
    f_1=0.
    f0=np.exp(-chi**2/2.)*(alpha_x/pi)**0.25
    fn=f0
    for i in range(1,n+1):
        fn=(chi*sq2*f0-sqrt(i-1)*f_1)/sqrt(i)
        f_1=f0
        f0=fn
    return fn

def harmonicfun_y(y,n):
    chi=sqrt(alpha_y)*y
    f_1=0.
    f0=np.exp(-chi**2/2.)*(alpha_y/pi)**0.25
    fn=f0
    for i in range(1,n+1):
        fn=(chi*sq2*f0-sqrt(i-1)*f_1)/sqrt(i)
        f_1=f0
        f0=fn
    return fn

#an array of harmonic functions to given order
def harmonicfun_x_all(x,nmax):
    fosc=np.zeros(nmax+1)
    chi=sqrt(alpha_x)*x
    f_1=0.
    f0=np.exp(-chi**2/2.)*(alpha_x/pi)**0.25
    fn=f0
    fosc[0]=fn
    for i in range(1,nmax+1):
        fn=(chi*sq2*f0-sqrt(i-1)*f_1)/sqrt(i)
        f_1=f0
        f0=fn
        fosc[i]=fn
    return fosc

def harmonicfun_y_all(y,nmax):
    fosc=np.zeros(nmax+1)
    chi=sqrt(alpha_y)*y
    f_1=0.
    f0=np.exp(-chi**2/2.)*(alpha_y/pi)**0.25
    fn=f0
    fosc[0]=fn
    for i in range(1,nmax+1):
        fn=(chi*sq2*f0-sqrt(i-1)*f_1)/sqrt(i)
        f_1=f0
        f0=fn
        fosc[i]=fn
    return fosc
    

def product_fun(x,y,i,j,ip,jp):
    vx_harm=.5*m_x*x*x*omega_x**2
    vy_harm=.5*m_y*y*y*omega_y**2
    vxy_correct=pot2D.potential_ho(x,y)
    vdiff=vxy_correct-vx_harm-vy_harm
    vg=vdiff*harmonicfun_x(x,i)*harmonicfun_x(x,ip)*harmonicfun_y(y,j)*\
    harmonicfun_y(x,jp)
    return vg
    
#Bound states in 2D potentials
import matplotlib.pyplot as plt
from mpl_toolkits  import mplot3d
import numpy as np
from numpy import linalg as LA
import scipy.integrate as integrate
import harmf,pot2D
#descretization in each dimension
ndivx=64
ndivy=64
#description of the basis
basis_size_x=5
basis_size_y=5
omega_x=1.
omega_y=1.
xmax,ymax=harmf.xy_dim(basis_size_x,basis_size_y)
#meshgrid
xarr=np.linspace(-xmax,xmax,num=ndivx)
yarr=np.linspace(-ymax,ymax,num=ndivy)
xx,yy=np.meshgrid(xarr,yarr)

#############################################################################
#now compute the Hamiltonian matrix
icmax=(basis_size_x+1)*(basis_size_y+1)
H_mat=np.zeros((icmax,icmax))
ic=-1
for i in range(0,basis_size_x+1):
    print(i)
    for j in range(0,basis_size_y+1):
        ic=ic+1
        icp=-1
        for ip in range(0,basis_size_x+1):
            for jp in range(0,basis_size_y+1):
                icp=icp+1
                if ic<=icp:                 
                    H_mat[ic,icp]=integrate.dblquad(harmf.product_fun,-xmax,\
                    xmax,-ymax,ymax,args=(i,j,ip,jp))[0]
                else:
                    H_mat[ic,icp]=H_mat[icp,ic]
        H_mat[ic,ic]=H_mat[ic,ic]+(i+.5)*omega_x+(j+.5)*omega_y
#find eigenvalues and eigenvectors of H_mat
eigen_E,eigen_V=LA.eig(H_mat)
idx=np.argsort(eigen_E)
eigen_E=eigen_E[idx]
eigen_V=eigen_V[:,idx]
#print the potential, ground state wavefunction and probability distribution
ifun=0
wave=np.zeros((ndivx,ndivy))
vxy=np.zeros((ndivx,ndivy))
for i in range(0,ndivx):
    for j in range(0,ndivy):
        fosc_x=harmf.harmonicfun_x_all(xx[i,j],basis_size_x)
        fosc_y=harmf.harmonicfun_y_all(yy[i,j],basis_size_y)
        s=0.
        ic=-1
        for nx in range(0,basis_size_x+1):
            for ny in range(0,basis_size_y+1):
                ic=ic+1
                s=s+eigen_V[ic,ifun]*fosc_x[nx]*fosc_y[ny]
        wave[i,j]=s
        vxy[i,j]=pot2D.potential_ho(xx[i,j],yy[i,j])
fig1=plt.figure(1,figsize=plt.figaspect(0.5))
fig1.tight_layout()
#plt.subplots_adjust(wspace=1)
ax=fig1.add_subplot(1,2,1,projection='3d')
ax.plot_surface(xx,yy,vxy,cmap='viridis')
plt.xlabel('x')
plt.ylabel('y')
plt.title('V(x,y)=2D-SH0+Gaussian')


ax=fig1.add_subplot(1,2,2,projection='3d')
ax.plot_surface(xx,yy,wave,cmap='plasma')
plt.xlabel('x')
plt.ylabel('y')
plt.title('\u03C6 (x,y)')
plt.savefig('2d_hog.png',bbox_inches = "tight")

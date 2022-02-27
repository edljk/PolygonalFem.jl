#########################################################################################
#########################################################################################
#####                Resolution of the linearized elasticity system                 #####
#####                      by the Virtual Element Method (VEM)                      #####
#########################################################################################
#########################################################################################

import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

#########################################################################################
#######             Boundary conditions in the case of a L-shaped domain          #######
#######             Input: c (n*2 array) Coordinates of point                     #######
#######             Output: g (n array) Corresponding values of bc                #######
#########################################################################################
def lshape_bc(c):
  
  EPS = 1e-4
  g = [0.0,0.0]
  isbc = 1 if ( c[1] > 1.0 - EPS ) else 0
    
  return (isbc,g)

#########################################################################################
#########################################################################################

#########################################################################################
#########################################################################################

#########################################################################################
#######             R.h.s. function in the case of a L-shaped domain              #######
#######             Input: c(table of 2 float) coordinates of point               #######
#######             Output: (float) value of r.h.s. at c                          #######
#########################################################################################
def lshape_rhs(c) :
  f = [0.0,0.0]
  
  if ( c[0] > 0.9 ) :
    f[1] = -1.0
   
  return f

#########################################################################################
#########################################################################################

#########################################################################################
#######        Calculate the volume, barycenter and diameter of an element        #######
#######             Input:   vc(table of N*2 float) coordinates of points         #######
#######             Outputs: vol (float) volume of elt                            #######
#######             Outputs: bar (list of size 2) barycenter of elt               #######
#######             Outputs: diam (float) diameter of elt                         #######
#########################################################################################
def characElt(vc) :

  n = vc.shape[0]

  # Size n array with entries x_iy_{i+1} - x_{i+1}y_i
  crossxiyi = vc[:,0] * np.roll(vc[:,1],-1) - np.roll(vc[:,0],-1) * vc[:,1]
  
  # vol = 0.5 * \sum (x_i y_{i+1} - x_{i+1}y_i )
  vol = 0.5 * np.abs(np.sum(crossxiyi))

  # bar_x = 1/(6*vol) * \sum (x_i + x_{i+1})(x_i y_{i+1} - x_{i+1}y_i )
  # bar_y = 1/(6*vol) * \sum (y_i + y_{i+1})(x_i y_{i+1} - x_{i+1}y_i )
  bar = 1.0/(6.0*vol) * np.sum( ( np.roll(vc,-1,axis=0) + vc ) * crossxiyi[:,None] , axis=0 )
    
  # Diameter
  diam = 0.0
  for i in range(n):
    for j in range(n):
      diam = max(diam,(vc[i,0]-vc[j,0])**2+ (vc[i,1]-vc[j,1])**2)
  diam = np.sqrt(diam)
  
  return (vol,bar,diam)
  
#########################################################################################
#########################################################################################



#########################################################################################
#######            Solve the Laplace equation                                    ########
#######            Inputs: smesh (string) .mat File of the mesh                  ########
#######                    rhs (function) Right-hand side function               ########
#######                    boundary_condition (function)                         ########
#######            Output: u (array with size 2 *np): solution at vertices       ########
#######                    u[2*ip+0] = u_x[ip] and u[2*ip+1] = u_y[ip]           ########
#########################################################################################
def vem(smesh,rhs,boundary_condition):

  # Material parameters
  lm = 0.5769
  mu = 0.3846

  # Load mesh
  mesh      = scipy.io.loadmat(smesh)
  vertices  = mesh['vertices']
  elements  = np.array([i[0].reshape(-1) - 1 for i in mesh['elements']],dtype=object)
  n         = vertices.shape[0]  # Total number of vertices
  ne        = elements.shape[0]  # Number of elements
  
  # Stiffness matrix, rhs and solution
  K = np.zeros((2*n,2*n))
  F = np.zeros(2*n)
  u = np.zeros(2*n)
  
  # Assembly of the local stiffness matrices and of the force vector
  for el in range(ne) :
  
    v     = elements[el] # Global indices of the vertices of el
    na    = v.shape[0]   # Number of vertices in el
    vc    = vertices[v]  # na * 2 array of vertices coordinates
    
    # Calculation of area, centroid and diameter of el
    (vol,bar,diam) = characElt(vc)
  
    # Local stiffness matrix
    # (3*3) Matrix D of elasticity over constant strain fields
    D = vol * np.array([
         [ 2.0*mu+lm ,     lm    ,   0.0  ],
         [     lm    , 2.0*mu+lm ,   0.0  ],
         [     0.0   ,   0.0     , 4.0*mu ]
        ])
    
    # Basis functions for R = [  (1,0)  ,  (0,1)  , (-(x_1-bar1,x_0-bar0)) ]
    # Basis functions for C = [ (x_0-bar0,0) , (0,x_1-bar1) , (x_1-bar1,x_0-bar0)  ]
    # (2n*3) Matrix WC of components of basis functions over constant strain fields
    WC = np.zeros((2*na,3))
    # (2n*3) Matrix WR of components of basis functions over rigid-body motions
    WR = np.zeros((2*na,3))
    # (2n*3) matrix NC: (NC)_{2i,l} = (c_k(x_i))_0, (NC)_{2i+1,l} = (c_k(x_i))_1
    NC = np.zeros((2*na,3))
    # (2n*3) matrix NR: (NR)_{2i,l} = (r_k(x_i))_0, (NR)_{2i+1,l} = (r_k(x_i))_1
    NR = np.zeros((2*na,3))
    
    # Loop over vertices to assemble WC and WR
    for i in range(na) :
      
      # Coordinates of vertices i-1, i, i+1
      cim = vc[(i-1) % na,:]
      ci  = vc[i,:]
      cip = vc[(i+1) % na,:]
  
      # Normal vector at vertex (multiplied by its length)
      nor = np.array([cip[1] - cim[1],-(cip[0] - cim[0])])
      
      q1 = 0.5*nor[0] / vol
      q2 = 0.5*nor[1] / vol
    
      WC[2*i,0]   = q1
      WC[2*i,1]   = 0.0
      WC[2*i,2]   = 0.5 * q2
      
      WC[2*i+1,0] = 0.0
      WC[2*i+1,1] = q2
      WC[2*i+1,2] = 0.5 * q1
      
      WR[2*i,0]   = 1.0 / na
      WR[2*i,1]   = 0.0
      WR[2*i,2]   = 0.5*q2
      
      WR[2*i+1,0] = 0.0
      WR[2*i+1,1] = 1.0 / na
      WR[2*i+1,2] = -0.5 * q1
      
      NC[2*i,0] = ci[0] - bar[0]
      NC[2*i,1] = 0.0
      NC[2*i,2] = ci[1] - bar[1]
        
      NC[2*i+1,0] = 0.0
      NC[2*i+1,1] = ci[1] - bar[1]
      NC[2*i+1,2] = ci[0] - bar[0]
      
      NR[2*i,0] = 1.0
      NR[2*i,1] = 0.0
      NR[2*i,2] = ci[1] - bar[1]
        
      NR[2*i+1,0] = 0.0
      NR[2*i+1,1] = 1.0
      NR[2*i+1,2] = -(ci[0] - bar[0])
      
    # Projection matrix over polynomial displacements
    PP = np.dot(WC,NC.T) + np.dot(WR,NR.T)
    
    # Projection matrix over higher-order functions
    QP = np.eye(2*na) - PP
    
    # Polynomial block in the local stiffness matrix
    PROJ = np.dot(WC,np.dot(D,WC.T))
    
    # Higher-order functions block in the local stiffness matrix
    alphaE = np.trace(D) / np.trace(np.dot(NC.T,NC))
    STAB = alphaE * np.dot(QP,QP.T)
    
    # Local stiffness matrix
    Ke = PROJ + STAB
    
    # Add local stiffness matrix to the global one
    for il in range(na) :
      ig = v[il]
      for jl in range(na) :
        jg = v[jl]
        K[2*ig,2*jg]     = K[2*ig,2*jg] + Ke[2*il,2*jl]
        K[2*ig,2*jg+1]   = K[2*ig,2*jg+1] + Ke[2*il,2*jl+1]
        K[2*ig+1,2*jg]   = K[2*ig+1,2*jg] + Ke[2*il+1,2*jl]
        K[2*ig+1,2*jg+1] = K[2*ig+1,2*jg+1] + Ke[2*il+1,2*jl+1]

    # Update force vector
    f = rhs(bar)
    for il in range(na) :
      ig = v[il]
      F[2*ig]    = F[2*ig] + vol / na * f[0]
      F[2*ig+1]  = F[2*ig+1] + vol / na * f[1]

  # Deal with boundary condition
  # Construct table of bc indices and bc values
  boundary = []
  bcval    = []
  for i in range(n) :
    c = vertices[i]
    (isbc,g) = boundary_condition(c)
    if isbc :
      boundary.extend([2*i,2*i+1])
      bcval.extend(g)
      
  boundary = np.array(boundary)
  bcval    = np.array(bcval)
  
  # Internal degrees of freedom
  idof  = np.array([ i for i in np.arange(2*n) if i not in boundary ]).ravel() # Table of internal dofs
  nidof = idof.shape[0]
  
  # Modification of the rhs
  F = F - np.dot(K[:,boundary],bcval)
  
  # Extraction of the sub-stiffness matrix K_II
  K_II = np.zeros((nidof,nidof))
  
  for il in range(nidof) :
    ig = idof[il]
    for jl in range(nidof) :
      jg = idof[jl]
      K_II[il,jl] = K[ig,jg]
  
  # Matrix inversion u_I = K_{II}^{—1} F_I
  u[idof]     = np.linalg.solve(K_II,F[idof])
  u[boundary] = bcval

  return u
    
#########################################################################################
#########################################################################################

#########################################################################################
#######            Plot function u defined at the vertices of mesh               ########
#######            Inputs: smesh (string) .mat File of the mesh                  ########
#######                    u (np array) Array containing sol.                    ########
#######                    save (Bool) Enable plot save                          ########
#######                    plot_name (string) Name of the plot                   ########
#########################################################################################
def plot_sol(smesh,u,save=False,plot_name=None) :

  # Load mesh
  mesh = scipy.io.loadmat(smesh)
  vertices = mesh['vertices']
  elements = np.array([i[0].reshape(-1) - 1 for i in mesh['elements']],dtype=object)
  
  ne = len(elements) # Number of elements

  x = vertices[:,0]
  y = vertices[:,1]
  vx = u[::2]
  vy = u[1::2]
  
  fig, axs = plt.subplots(1,2,figsize=(12,5),constrained_layout=False)
  
  xi = np.linspace(min(x) - 0.01, max(x) + 0.001, 100)
  yi = np.linspace(min(y) - 0.01, max(y) + 0.001, 100)
  
  # Interpolate values of u_x and u_y at the grid nodes
  zxi = griddata((x, y), vx, (xi[None, :], yi[:, None]), method='linear')
  zyi = griddata((x, y), vy, (xi[None, :], yi[:, None]), method='linear')

  # Set title and axes
  axs[0].set_title("Approximate displacement (u_x)")
  axs[1].set_title("Approximate displacement (u_y)")
  axs[0].set_xlabel('x')
  axs[0].set_ylabel('y')
  axs[1].set_xlabel('x')
  axs[1].set_ylabel('y')
  
  # Plot elements
  for i in range(ne) :
    npe = len(elements[i])
    for j in range(npe) :
      ip0 = elements[i][j]
      ip1 = elements[i][(j+1)%npe]
      xx = [x[ip0],x[ip1]]
      yy = [y[ip0],y[ip1]]
      axs[0].plot(xx,yy,"k",linewidth=0.5)
      axs[1].plot(xx,yy,"k",linewidth=0.5)
  
  # Plot solution
  imx = axs[0].pcolormesh(xi, yi, zxi)
  imy = axs[1].pcolormesh(xi, yi, zyi)

  divider = make_axes_locatable(axs[0])
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(imx, cax=cax)

  divider = make_axes_locatable(axs[1])
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(imy, cax=cax)

  # Save plot if need be
  if save and plot_name is not None:
      plt.savefig(plot_name)
  elif save and plot_name is None:
      plt.savefig("sol.png")
  
  fig.subplots_adjust(wspace=0.5)
  plt.show()

#########################################################################################
#########################################################################################

#########################################################################################
################                       Main function                     ################
#########################################################################################
def main():
  import argparse
  from argparse import RawTextHelpFormatter
  
  parser = argparse.ArgumentParser(description='Solver for the 2d Poisson equation on polygonal elements', formatter_class = RawTextHelpFormatter)
  parser.add_argument("-i", help = "Path to input mesh", type=str)
  parser.add_argument("-o", help = "Path to output solution", type=str, default="./sol.npy")
  parser.add_argument("-d", help = "Shape of domain:\n- s: Square domain\n- l: L-shaped domain",type=str)
  parser.add_argument("--save_plot", help="Enable save plot", action="store_true")
  parser.add_argument("--title", help="Title of the plot", type=str, default="./plot.png")
  
  args = parser.parse_args()
  smesh = args.i
  
  # Vector for solution
  u = None
  
  # Resolution of the system
  u = vem(smesh,lshape_rhs,lshape_bc)

  # Save file
  # np.save(args.o,u)
  
  # Plot mesh and solution
  plot_sol(smesh,u,args.save_plot,args.title)
  
#########################################################################################
#########################################################################################

# Call main when file is executed as a script
if __name__ == '__main__':
  main()

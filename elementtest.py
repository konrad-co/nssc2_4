from stiffnessmatrixgeneration import *
import matplotlib.pyplot as plt
from importdata import *


N=10

mesh = generateMesh(L,hz,N)


Dt=BoundaryCondition('top', 'Dirichlet' ,bc_top)
Nb=BoundaryCondition('bottom', 'Neumann', bc_bottom)
Nl=BoundaryCondition('left', 'Neumann', 0)
Nr=BoundaryCondition('right', 'Neumann', 0)

sol=solveBVP(mesh, [Dt,Nb,Nl,Nr], k)


#calculate gradients and fluxes and store them for each triangle
calcLocalGradsAndFluxes(mesh,sol,k)


'''
for triangle in mesh.elements:
	print(np.sum(triangle.tempgrad))'''


from stiffnessmatrixgeneration import *
import matplotlib.pyplot as plt
from importdata import *


N=10
#print(elements_modify)

mesh = generateMesh(L,hz,N)


Dt=BoundaryCondition('top', 'Dirichlet' ,bc_top)
Nb=BoundaryCondition('bottom', 'Neumann', bc_bottom/k)
Nl=BoundaryCondition('left', 'Neumann', 0)
Nr=BoundaryCondition('right', 'Neumann', 0)


sol=solveBVP(mesh, [Dt,Nb,Nl,Nr], k, 0, elements_modify, c)


#calculate gradients and fluxes and store them for each triangle
calcLocalGradsAndFluxes(mesh,sol,k,4,elements_modify,c)


'''
for triangle in mesh.elements:
	print(np.sum(triangle.tempgrad))'''


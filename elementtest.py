from stiffnessmatrixgeneration import *
#from importdata import *

L = 1
k = 1
hz = 0.001
bc_top=100
bc_bottom = 100.
N=4

mesh = generateMesh(L,hz,N)

rhs=np.zeros(mesh.number_elements)

Dt=BoundaryCondition('top', 'Dirichlet' ,bc_top)
Nb=BoundaryCondition('bottom', 'Neumann', bc_bottom)
Nl=BoundaryCondition('left', 'Neumann', 0)
Nr=BoundaryCondition('right', 'Neumann', 0)

sol=solveBVP(mesh, [Dt,Nb,Nl,Nr], k)

#print(sol)

calcLocalGradsAndFluxes(mesh,sol,k)

#for triangle in mesh.elements:
	#print(np.sum(triangle.flux))


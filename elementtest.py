from stiffnessmatrixgeneration import *
from importdata import *

N=10

mesh = generateMesh(L,hz,N)

Dt=BoundaryCondition('top', 'Dirichlet' ,bc_top)
Nb=BoundaryCondition('bottom', 'Neumann', bc_bottom)
Nl=BoundaryCondition('left', 'Neumann', 0)
Nr=BoundaryCondition('right', 'Neumann', 0)

sol=solveBVP(mesh, [Dt,Nb,Nl,Nr], k)

print(sol)

calcLocalGradsAndFluxes(mesh,sol,k)



for triangle in mesh.elements:
	print(np.sum(triangle.flux))


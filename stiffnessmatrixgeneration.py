import numpy as np
from meshgeneration import *


def ElementStiffnessMatrix(Triangle, k,N, L):
	'''Input: 	Triangle	element of the class Triangle
				k			conductivity associated with the triangle
				N 			number of gridpoints per direction
				L			length of the domain
				hz			element thickness
	'''
	vertex1=nodeToCoordinate(Triangle.nodes[0],N,L)
	vertex2=nodeToCoordinate(Triangle.nodes[1],N,L)
	vertex3=nodeToCoordinate(Triangle.nodes[2],N,L)


	b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])
	c=np.array([vertex3[0]-vertex2[0], vertex1[0]-vertex3[0], vertex2[0]-vertex1[0]])

	He=k/(4*Triangle.area)*Triangle.hz*(np.outer(b,b)+np.outer(c,c))
	return He

def assembleGlobalStiffnessMatrix(mesh,k):
	'''Input:	mesh 		underlying mesh
				k			conductivity
	'''
	H=np.zeros((mesh.N*mesh.N, mesh.N*mesh.N))
	'''for element in mesh.elements:
		He=ElementStiffnessMatrix(element, k, mesh.N, mesh.L)
		H[np.ix_([element.nodes],element.nodes)] += He'''
	for element in mesh.elements:
		He=ElementStiffnessMatrix(element, k, mesh.N, mesh.L)
		for i in range(np.shape(He)[0]):
			for j in range(np.shape(He)[1]):
				H[element.nodes[i],element.nodes[j]]=He[i,j]
		
	return H



def solveBVP(mesh, boundaries,k):
	'''Input: 		mesh 		# element of class Mesh
					rhs			# vector containing the right hand side of the system
					boundaries 	# list of elements of the class BoundaryCondition
	'''

	#assemble stiffnessmatrix
	H=assembleGlobalStiffnessMatrix(mesh,k)
	#print('H=',H)

	#get Dirichlet boundary nodes
	bottom, right, top, left = mesh.getBoundary()
	dirichletnodes=[]
	T=0
	for element in boundaries:
		if element.type == 'Dirichlet':
			if element.location == 'left':
				dirichletnodes.extend(left)
				T=element.value
			if element.location == 'right':
				dirichletnodes.extend(right)
				T=element.value
			if element.location == 'bottom':
				dirichletnodes.extend(bottom)
				T=element.value
			if element.location == 'top':
				dirichletnodes.extend(top)
				T=element.value
	dirichletnodes=np.unique(np.array(dirichletnodes))

	#get Neumann nodes
	neumannnodes=[]
	q=0
	for element in boundaries:
		if element.type=='Neumann' and element.value != 0:
			if element.location == 'left':
				neumannnodes.extend(left)
				q=element.value
			if element.location == 'right':
				neumannnodes.extend(right)
				q=element.value
			if element.location == 'bottom':
				neumannnodes.extend(bottom)
				q=element.value
			if element.location == 'top':
				neumannnodes.extend(top)
				q=element.value


	rhs=np.zeros(mesh.N*mesh.N)
	for node in neumannnodes:
		rhs[node]-=nodalForce(mesh,q,k)

	#extract free nodes
	freenodes=list(set(range(mesh.N*mesh.N)).difference(dirichletnodes))
	

	A=H[np.ix_(freenodes,freenodes)]
	print(A)

	
	#build RHS including Neumann BC
	#print(rhs[freenodes])
	rhs[freenodes]=rhs[freenodes] - H[np.ix_(freenodes,dirichletnodes)] @ (np.ones(len(dirichletnodes))*T)
	#print(np.shape(rhs_subsystem))
	print('rhs=',rhs)
	#print(H[np.ix_(freenodes,dirichletnodes)] @ (np.ones(len(dirichletnodes))*T))
	

	#solve the system

	sol = np.ones(mesh.N*mesh.N)*T
	#sol[freenodes] = np.linalg.solve(A, rhs_subsystem)
	sol[freenodes]=np.linalg.solve(A, rhs[freenodes])
	#print(sol[freenodes])

	#compute the remaining reaction forces... why?
	Prf=H[np.ix_(dirichletnodes,freenodes)] @ sol[freenodes] + H[np.ix_(dirichletnodes,dirichletnodes)] @ sol[dirichletnodes]

	return sol

def nodalForce(mesh,q,k):
	return mesh.elements[0].hz*mesh.L/(mesh.N-1)*q/2

def calcLocalGradsAndFluxes(mesh,sol,k):
	for triangle in mesh.elements:
		vertex1=nodeToCoordinate(triangle.nodes[0],mesh.N,mesh.L)
		vertex2=nodeToCoordinate(triangle.nodes[1],mesh.N,mesh.L)
		vertex3=nodeToCoordinate(triangle.nodes[2],mesh.N,mesh.L)

		b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])
		c=np.array([vertex3[0]-vertex2[0], vertex1[0]-vertex3[0], vertex2[0]-vertex1[0]])

		grad=1/(2*triangle.area)*np.array([b,c]) @ sol[triangle.nodes]
		triangle.tempgrad=grad
		triangle.flux=-k*grad






	








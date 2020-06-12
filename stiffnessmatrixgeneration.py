import numpy as np
from meshgeneration import *


def ElementStiffnessMatrix(Triangle, k,N, L):
	'''Input: 	Triangle	element of the class Triangle
				k			conductivity associated with the triangle
				N 			number of gridpoints per direction
				L			length of the domain
				hz			element thickness
	'''
	vertex1, vertex2, vertex3=Triangle.nodes
	vertex1=nodeToCoordinate(vertex1,N,L)
	vertex2=nodeToCoordinate(vertex2,N,L)
	vertex3=nodeToCoordinate(vertex3,N,L)


	a=np.array([vertex2[0]*vertex3[1]-vertex3[0]*vertex2[1], vertex3[0]*vertex1[1]-vertex1[0]*vertex3[1], vertex1[0]*vertex2[1]-vertex2[0]*vertex1[0]])
	b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])
	c=np.array([vertex3[0]-vertex2[0], vertex1[0]-vertex3[0], vertex2[0]-vertex1[0]])

	area=vertex1[0]*b[0]+vertex2[0]*b[1]+vertex3[0]*b[2]

	#He=k/(4*area*Triangle.hz)*np.array([[b[0]*b[0]+c[0]*c[0], b[0]*b[1]+c[0]*c[1], b[0]*b[2]+c[0]*c[2]], [b[1]*b[0]+c[1]*c[0], b[1]*b[1]+c[1]*c[1], b[1]*b[2]+c[1]*c[2]], [b[2]*b[0]+c[2]*c[0], b[2]*b[1]+c[2]*c[1], b[2]*b[2]+c[2]*c[2]]])
	He=k/(4*area*Triangle.hz)*(np.outer(b,b)+np.outer(c,c))
	return He

def assembleGlobalStiffnessMatrix(mesh,k):
	'''Input:	mesh 		underlying mesh
				k			conductivity
	'''
	H=np.zeros((mesh.number_elements, mesh.number_elements))
	for element in mesh.elements:
		He=ElementStiffnessMatrix(element, k, mesh.N, mesh.L)
		H[np.ix_(element.nodes,element.nodes)] += He
	return H

def solveBVP(mesh, boundaries,k):
	'''Input: 		mesh 		# element of class Mesh
					rhs			# vector containing the right hand side of the system
					boundaries 	# list of elements of the class BoundaryCondition
	'''

	#assemble stiffnessmatrix
	H=np.zeros((mesh.N*mesh.N, mesh.N*mesh.N))
	for element in mesh.elements:
		He=ElementStiffnessMatrix(element, k, mesh.N, mesh.L)
		#print(He)
		H[np.ix_(element.nodes,element.nodes)] += He
		#print(np.ix_(element.nodes,element.nodes))

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

	#rhs=np.ones(mesh.N*mesh.N)
	rhs=np.zeros(mesh.N*mesh.N)
	for node in neumannnodes:
		rhs[node]-=nodalForce(mesh,q)

	#extract free nodes
	freenodes=list(set(range(mesh.N*mesh.N)).difference(dirichletnodes))
	#print(freenodes)

	A=H[np.ix_(freenodes,freenodes)]

	#build RHS including Neumann BC
	rhs[freenodes]=rhs[freenodes] - H[np.ix_(freenodes,dirichletnodes)] @ (np.ones(len(dirichletnodes))*T)
	#print(rhs)
	#print(H[np.ix_(freenodes,dirichletnodes)] @ (np.ones(len(dirichletnodes))*T))

	#solve the system
	sol = np.ones(mesh.N*mesh.N)*T
	sol[freenodes] = np.linalg.solve(A, rhs[freenodes])

	#compute the remaining reaction forces... why?
	Prf=H[np.ix_(dirichletnodes,freenodes)] @ sol[freenodes] + H[np.ix_(dirichletnodes,dirichletnodes)] @ sol[dirichletnodes]

	return sol

def nodalForce(mesh,q):
	return mesh.elements[0].hz*mesh.L/(mesh.N-1)*q/2

def calcLocalGradsAndFluxes(mesh,sol,k):
	for triangle in mesh.elements:
		vertex1, vertex2, vertex3=triangle.nodes
		vertex1=nodeToCoordinate(triangle.nodes[0],mesh.N,mesh.L)
		vertex2=nodeToCoordinate(triangle.nodes[1],mesh.N,mesh.L)
		vertex3=nodeToCoordinate(triangle.nodes[2],mesh.N,mesh.L)

		b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])
		c=np.array([vertex3[0]-vertex2[0], vertex1[0]-vertex3[0], vertex2[0]-vertex1[0]])

		area=vertex1[0]*b[0]+vertex2[0]*b[1]+vertex3[0]*b[2]

		grad=1/(2*area)*np.array([b,c]) @ sol[triangle.nodes]
		triangle.tempgrad=grad
		triangle.flux=-k*grad


	








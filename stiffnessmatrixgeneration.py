import numpy as np
from meshgeneration import *
from print_HTP import *


def ElementStiffnessMatrix(Triangle, k,N, L, V=0):
	'''Input: 	Triangle	element of the class Triangle
				k			conductivity associated with the triangle
				N 			number of gridpoints per direction
				L			length of the domain
				hz			element thickness
	'''
	vertex1=nodeToCoordinate(Triangle.nodes[0],N,L, V)
	vertex2=nodeToCoordinate(Triangle.nodes[1],N,L, V)
	vertex3=nodeToCoordinate(Triangle.nodes[2],N,L, V)


	b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])
	c=np.array([vertex3[0]-vertex2[0], vertex1[0]-vertex3[0], vertex2[0]-vertex1[0]])

	He=k/(4*Triangle.area)*Triangle.hz*(np.outer(b,b)+np.outer(c,c))
	return He

def assembleGlobalStiffnessMatrix(mesh,k, V=0, elements_modify = [], c = 1):
	'''Input:	mesh 		underlying mesh
				k			conductivity
	'''

	# Test variations of coordinates
	#print("V=0", nodeToCoordinate(99, mesh.N, mesh.L, V=0))
	#print("V=3", nodeToCoordinate(99, mesh.N, mesh.L, V=3))
	#print("V=0", nodeToCoordinate(90, mesh.N, mesh.L, V=0))
	#print("V=3", nodeToCoordinate(90, mesh.N, mesh.L, V=3))

	H=np.zeros((mesh.N*mesh.N, mesh.N*mesh.N))
	for element in mesh.elements:
		if V == 4 and element.id in elements_modify:
			He = ElementStiffnessMatrix(element, k*c, mesh.N, mesh.L, V)
		else:
			He = ElementStiffnessMatrix(element, k, mesh.N, mesh.L, V)
		H[np.ix_(element.nodes,element.nodes)] += He		
	return H




def solveBVP(mesh, boundaries,k, V=0, elements_modify = [], c = 1):
	'''Input: 		mesh 		# element of class Mesh
					rhs			# vector containing the right hand side of the system
					boundaries 	# list of elements of the class BoundaryCondition
	'''

	#assemble stiffnessmatrix
	H=assembleGlobalStiffnessMatrix(mesh,k, V, elements_modify, c)

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

	#get Neumann boundarynodes
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


	#extract free nodes
	freenodes=list(set(range(mesh.N*mesh.N)).difference(dirichletnodes))

	
	#build RHS including Neumann BC
	rhs=np.zeros(mesh.N*mesh.N)
	for node in neumannnodes:
		if V == 4 & node in elements_modify:
			rhs[node]-=nodalForce(mesh,q,k*c)
		else:
			rhs[node]-=nodalForce(mesh,q,k)
	

	#solve the system
	sol = np.ones(mesh.N*mesh.N)*T
	sol[freenodes]=np.linalg.solve(H[np.ix_(freenodes,freenodes)], rhs[freenodes]- H[np.ix_(freenodes,dirichletnodes)] @ (np.ones(len(dirichletnodes))*T))


	#compute the remaining reaction forces and store them in rhs
	rhs[dirichletnodes]=H[np.ix_(dirichletnodes,freenodes)] @ sol[freenodes] + H[np.ix_(dirichletnodes,dirichletnodes)] @ sol[dirichletnodes]

	#print to output file
	print_HTP(H,sol.reshape(mesh.N,mesh.N),rhs.reshape(mesh.N,mesh.N))

	return sol

def nodalForce(mesh,q,k):
	return mesh.elements[0].hz*mesh.L/(mesh.N-1)*q/2

def calcLocalGradsAndFluxes(mesh,sol,k,V=0, elements_modif=[], c=1):
	for triangle in mesh.elements:
		vertex1=nodeToCoordinate(triangle.nodes[0],mesh.N,mesh.L)
		vertex2=nodeToCoordinate(triangle.nodes[1],mesh.N,mesh.L)
		vertex3=nodeToCoordinate(triangle.nodes[2],mesh.N,mesh.L)

		b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])
		cvec=np.array([vertex3[0]-vertex2[0], vertex1[0]-vertex3[0], vertex2[0]-vertex1[0]])



		grad=1/(2*triangle.area)*np.array([b,cvec]) @ sol[triangle.nodes]
		triangle.tempgrad=grad
		if V==4 and element in elements_modify:
			triangle.flux=-k*c*grad
		else:
			triangle.flux=-k*grad






	








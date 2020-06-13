#code is based on the solution of an exerise for the course "Numerik partieller Differentialgleichungen: Stationäre Probleme" by Prof. Feischl
import numpy as np


class Triangle:
	def __init__(self, id, nodes, bndryedges, hz, tempgrad=np.zeros(2), flux=np.zeros(2), area=0):
		self.id=id 					# unique identification number of the element according to Figure 1 of the handout
		self.nodes=nodes			# indices of the vertíces
		self.bndryedges=bndryedges 	# array which specifies if the edges [1,2], [2,3], [3,1] are at the bottom/right/top/left boundary or inside the domain
		self.hz=hz					# element thickness
		self.tempgrad=tempgrad		# gradient of the temperature
		self.flux=flux				# flux of the temperature
		self.area=area

class Mesh:
	def __init__(self, L, elements, N):
		self.L=L 								# length of the domain
		self.elements=elements 					# list of triangles
		self.number_elements=len(elements)
		self.N=N 
		calcArea(self)								# number of gridpoints in each direction

	def getBoundary(self):
		'''returns 4 lists containing the nodes of each boundary'''
		leftboundary=[]
		rightboundary=[]
		topboundary=[]
		bottomboundary=[]

		for element in self.elements:
			if len(set(element.bndryedges).intersection({'left', 'right','top','bottom'})) != 0:

				if element.bndryedges[0] == 'left':
					leftboundary.extend([element.nodes[0], element.nodes[1]])
				if element.bndryedges[1] == 'left':
					leftboundary.extend([element.nodes[1], element.nodes[2]])
				if element.bndryedges[2] == 'left':
					leftboundary.extend([element.nodes[2], element.nodes[0]])

				if element.bndryedges[0] == 'right':
					rightboundary.extend([element.nodes[0], element.nodes[1]])
				if element.bndryedges[1] == 'right':
					rightboundary.extend([element.nodes[1], element.nodes[2]])
				if element.bndryedges[2] == 'right':
					rightboundary.extend([element.nodes[2], element.nodes[0]])

				if element.bndryedges[0] == 'top':
					topboundary.extend([element.nodes[0], element.nodes[1]])
				if element.bndryedges[1] == 'top':
					topboundary.extend([element.nodes[1], element.nodes[2]])
				if element.bndryedges[2] == 'top':
					topboundary.extend([element.nodes[2], element.nodes[0]])

				if element.bndryedges[0] == 'bottom':
					bottomboundary.extend([element.nodes[0], element.nodes[1]])
				if element.bndryedges[1] == 'bottom':
					bottomboundary.extend([element.nodes[1], element.nodes[2]])
				if element.bndryedges[2] == 'bottom':
					bottomboundary.extend([element.nodes[2], element.nodes[0]])

		return bottomboundary,rightboundary,topboundary,leftboundary
		#return np.unique(np.array(bottomboundary)), np.unique(np.array(rightboundary)), np.unique(np.array(topboundary)), np.unique(np.array(leftboundary))


class BoundaryCondition:
	def __init__(self, location, type, value):
		self.location=location
		self.type=type		
		self.value=value




def nodeToCoordinate(node,N,L):
	j=(node%N)
	i=(node-j)/N
	return np.array([i*L/(N-1),j*L/(N-1)])




def generateMesh(L, hz, N=10):
	#N number of nodes in one direction, N=10 by default

	#generate Triangles
	triangle_id=0
	elements=[]
	node_number=0		#node_number=i*N+j


	#row at the bottom
	#lower left corner
	elements.extend([Triangle(triangle_id,[node_number, node_number+1, node_number+N],['bottom','inner','left'], hz), Triangle(triangle_id+1,[node_number+1, node_number+1+N, node_number+N],['inner','inner','inner'], hz)])
	triangle_id+=2

	for j in range(1,N-2):
		meshelement_lower=Triangle(triangle_id,[j, j+1, j+N],['bottom','inner','inner'], hz)
		triangle_id+=1
		meshelement_upper=Triangle(triangle_id,[j+1, j+1+N, j+N], ['inner','inner','inner'], hz)
		triangle_id+=1
		elements.extend([meshelement_lower, meshelement_upper])
	
	#lower right corner
	node_number=N-2
	meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['bottom','inner','inner'], hz)
	triangle_id += 1
	meshelement_upper = Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N],['right','inner','inner'], hz)
	triangle_id += 1
	elements.extend([meshelement_lower, meshelement_upper])


	#inner nodes
	for i in range(1,N-2):
		#at the left boundary
		node_number=i*N
		meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['inner','inner','left'], hz)
		triangle_id+=1
		meshelement_upper=Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N],['inner','inner','inner'], hz)
		triangle_id+=1
		elements.extend([meshelement_lower, meshelement_upper])

		for j in range(1,N-2):
			node_number=i*N+j
			meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['inner','inner','inner'], hz)
			triangle_id+=1
			meshelement_upper=Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N], ['inner','inner','inner'], hz)
			triangle_id+=1
			elements.extend([meshelement_lower, meshelement_upper])

		#at the right boundary	
		node_number=i*N+(N-2)
		meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['inner','inner','inner'], hz)
		triangle_id += 1
		meshelement_upper = Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N],['right','inner','inner'], hz)
		triangle_id += 1
		elements.extend([meshelement_lower, meshelement_upper])


	#row at the top
	i=N-2
	#upper left corner
	node_number=i*N
	meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['inner','inner','left'], hz)
	triangle_id+=1
	meshelement_upper=Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N],['inner','top','inner'], hz)
	triangle_id+=1
	elements.extend([meshelement_lower, meshelement_upper])

	for j in range(1,N-2):
		node_number=i*N+j
		meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['inner','inner','inner'], hz)
		triangle_id+=1
		meshelement_upper=Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N],['inner','top','inner'], hz)
		triangle_id+=1
		elements.extend([meshelement_lower, meshelement_upper])

	#upper right corner
	node_number=i*N+N-2
	meshelement_lower=Triangle(triangle_id,[node_number, node_number+1, node_number+N],['inner','inner','inner'], hz)
	triangle_id+=1
	meshelement_upper=Triangle(triangle_id,[node_number+1, node_number+1+N, node_number+N],['right','top','inner'], hz)
	triangle_id+=1
	elements.extend([meshelement_lower, meshelement_upper])

	return Mesh(L,elements,N)

def calcArea(mesh):
	for triangle in mesh.elements:
		vertex1=nodeToCoordinate(triangle.nodes[0],mesh.N,mesh.L)
		vertex2=nodeToCoordinate(triangle.nodes[1],mesh.N,mesh.L)
		vertex3=nodeToCoordinate(triangle.nodes[2],mesh.N,mesh.L)

		b=np.array([vertex2[1]-vertex3[1], vertex3[1]-vertex1[1], vertex1[1]-vertex2[1]])

		triangle.area=abs(vertex1[0]*b[0]+vertex2[0]*b[1]+vertex3[0]*b[2])/2





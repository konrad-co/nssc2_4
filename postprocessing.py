from elementtest import *
import matplotlib.pyplot as plt
import matplotlib.lines as lines

#useful data
coordinates=nodeToCoordinate(np.array(range(N*N)).reshape((N,N)),N,L)
#edges=np.ix_(coordinates[0],coordinates[0])
#print(edges)
print(np.shape(coordinates))
print(coordinates)

'''
fig,ax = plt.subplots()
#ax.set_xticklabels(np.linspace(0,L,N))
#ax.set_yticklabels(np.linspace(0,L,N))
ax.set_xticks(np.linspace(0,L,N), minor=True)
ax.set_yticks(np.linspace(0,L,N), minor=True)
ax.set(xlim=(0,L), ylim=(0,L))
ax.set_aspect('equal')
plt.grid(which='minor')

#plt.hlines(np.linspace(0,L,N), 0, L,linewidths=0.4)
#plt.vlines(np.linspace(0,L,N), 0, L,linewidths=0.4)
plt.scatter(coordinates[0],coordinates[1], s=4, color='black')

plt.show()
'''

#plotting the temperature field
fig,ax = plt.subplots()
ax.set_xticks(np.linspace(0,L,N))
ax.set_yticks(np.linspace(0,L,N))
ax.set(xlim=(0,L), ylim=(0,L))
ax.set_aspect('equal')
#plt.grid(xdata=np.linspace(0,L,N, endpoint=True),ydata=np.linspace(0,L,N, endpoint=True))
plt.pcolormesh(coordinates[0],coordinates[1],sol.reshape((N,N),order='F'))
fig.colorbar(plt.pcolormesh(coordinates[0],coordinates[1],sol.reshape((N,N),order='F')))
plt.scatter(coordinates[0],coordinates[1], s=6, color='grey')
plt.xlabel("x in m")
plt.ylabel('y in m')
plt.title("Temperature field in K")
plt.savefig('temperaturefield.png')




#calculating the element centroids
centroids=[]
for element in mesh.elements:
	centroids.append((nodeToCoordinate(element.nodes[0],N,L)+nodeToCoordinate(element.nodes[1],N,L)+nodeToCoordinate(element.nodes[2],N,L))/3)

fig,ax = plt.subplots()
ax.set_xticks(np.linspace(0,L,N))
ax.set_yticks(np.linspace(0,L,N))
ax.set(xlim=(0,L), ylim=(0,L))
ax.set_aspect('equal')
for element in mesh.elements:
	#plt.quiver([centroids[element.id][0],centroids[element.id][1]],element.tempgrad[0],element.tempgrad[1])
plt.xlabel("x in m")
plt.ylabel('y in m')
plt.title("Gradient field")
plt.savefig('gradientfield.png')
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from src.OBJIO import loadOBJ
from src.fitting import getBounds, fitSphere, distance

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please provide an OBJ file as an argument'
		sys.exit(0)

	radiusNominal = 150

	vertices, faces, normals = loadOBJ(sys.argv[1])

	indicesLeft = np.where(vertices.T[0] < 0)
	indicesRight = np.where(vertices.T[0] > 0)

	verticesLeft = vertices[indicesLeft]
	verticesRight = vertices[indicesRight]

	boundsLeft = getBounds(verticesLeft.T)
	print boundsLeft
	boundsRight = getBounds(verticesRight.T)
	print boundsRight
	
	p0Left = [boundsLeft[0][0],boundsLeft[1][0],boundsLeft[2][0],radiusNominal]
	p0Right = [boundsRight[0][0],boundsRight[1][0],boundsRight[2][0],radiusNominal]
	centerPointLeft, radiusLeft = fitSphere(verticesLeft, p0Left, radiusNominal, boundsLeft)
	centerPointRight, radiusRight = fitSphere(verticesRight, p0Right, radiusNominal, boundsRight)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.scatter(verticesLeft.T[0], verticesLeft.T[1], verticesLeft.T[2], color='r')
	ax.scatter(centerPointLeft[0], centerPointLeft[1], centerPointLeft[2], color='yellow')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')

	ax.scatter(verticesRight.T[0], verticesRight.T[1], verticesRight.T[2], color='b')
	ax.scatter(centerPointRight[0], centerPointRight[1], centerPointRight[2], color='yellow')
	ax.set_xlabel('x')

	plt.show()
	
	print 'Distance: ' + str(distance(centerPointLeft, centerPointRight))

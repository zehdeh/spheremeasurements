import numpy as np

def getVoxelCenter(gridSize, gridScale, i, j, k):
	x = (i*gridScale[0]) - (gridSize[0]*gridScale[0]/2) + gridScale[0]/2
	y = (j*gridScale[1]) - (gridSize[1]*gridScale[1]/2) + gridScale[1]/2
	z = (k*gridScale[2]) - (gridSize[2]*gridScale[2]/2) + gridScale[2]/2

	return x, y, z

def getGridIndices(gridSize, gridScale, vertex):
	x = vertex[0] + gridSize[0]*gridScale[0]/2
	y = vertex[1] + gridSize[1]*gridScale[1]/2
	z = vertex[2] + gridSize[2]*gridScale[2]/2

	k = int(x / gridScale[0])
	j = int(y / gridScale[1])
	i = int(z / gridScale[2])

	return i,j,k

def isInGrid(gridSize, i, j, k):
	return i > 0 and i < gridSize[0] and j > 0 and j < gridSize[1] and k > 0 and k < gridSize[2]

def generateErrorGrid(gridSize, gridScale, spheres, verbose=False):
	return generateGrid(gridSize, gridScale, spheres, lambda(x): np.abs(x.fittingError()), verbose)

def generateCurvatureGrid(gridSize, gridScale, spheres, verbose=False):
	return generateGrid(gridSize, gridScale, spheres, lambda(x): x.curvature, verbose)
	
def generateGrid(gridSize, gridScale, spheres, evaluationFunction, verbose=False):
	centerPoints = []
	totalErrorMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
	nMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.int32)

	maxError = -np.inf

	for sphere in spheres:
		errors = evaluationFunction(sphere)
		errorMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
		for vertex,error in zip(sphere.vertices, errors):
			i,j,k = getGridIndices(gridSize, gridScale, vertex)

			if error > maxError:
				maxError = error
				sphereWithMaxError = sphere.fileName

			if isInGrid(gridSize,i,j,k):
				errorMatrix[i,j,k] += error
				nMatrix[i,j,k] += 1
			else:
				raise RuntimeError('There are points outside your grid. Is it too small?')
		totalErrorMatrix += errorMatrix

	if verbose:
		print sphereWithMaxError

	# Average over all samples
	totalErrorMatrix[np.nonzero(nMatrix)] = totalErrorMatrix[np.nonzero(nMatrix)] / nMatrix[np.nonzero(nMatrix)]

	# Normalize to 1.0
	totalErrorMatrix[np.nonzero(nMatrix)] = totalErrorMatrix[np.nonzero(nMatrix)] / totalErrorMatrix.max()
	
	return totalErrorMatrix, nMatrix

def generateVectorField(gridSize, gridScale, spheres, nMatrix, errorMatrix=None):
	totalVectorField = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=(np.float,3))

	if errorMatrix is None:
		errorMatrix = np.ones((gridSize[0], gridSize[1], gridSize[2]))
	for sphere in spheres:
		vectorField = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=(np.float,3))
		for vertex,normal in zip(sphere.vertices, sphere.normals):
			mag = np.linalg.norm(normal)
			if mag == 0:
				continue
			i,j,k = getGridIndices(gridSize, gridScale, vertex)

			error = errorMatrix[i,j,k]
			if isInGrid(gridSize,i,j,k):
				normal = normal/mag
				vectorField[i,j,k] += normal*error
		totalVectorField += vectorField
	totalVectorField[np.nonzero(nMatrix)] = totalVectorField[np.nonzero(nMatrix)] / nMatrix[np.nonzero(nMatrix)][...,None]
	return totalVectorField



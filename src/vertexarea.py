import numpy as np
import math
from timeit import default_timer as timer

def getFaceAngles(face, vertices):
	v = vertices[face]
	a = np.linalg.norm(v[0] - v[1])
	b = np.linalg.norm(v[1] - v[2])
	c = np.linalg.norm(v[2] - v[0])

	B = math.acos((a**2 + c**2 - b**2) / (2*a*c))
	C = math.asin((c*math.sin(B))/b)
	A = np.pi - B - C
	return A,B,C

def getFaceArea(face, vertices):
	v = vertices[face]
	a = np.linalg.norm(v[0] - v[1])
	b = np.linalg.norm(v[1] - v[2])
	c = np.linalg.norm(v[2] - v[0])

	t = [a,b,c]

	s = np.sum(t)/2
	return np.sqrt(s*np.prod(s-t))

def getVertexAreas(faces, vertices):
	faceAreas = np.array([getFaceArea(f, vertices) for f in faces])/3

	vertexAreas = np.array([np.sum(faceAreas[np.argwhere(np.any(np.reshape(i == faces.ravel(), (-1,3)), axis=1))]) for i,v in enumerate(vertices)])

	return vertexAreas

import numpy as np
from scipy.optimize import leastsq, least_squares

def distance(p1,p2):
	return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2)

def fittingErrorSphere(center, vertices):
	x0,y0,z0,R = center
	x,y,z = vertices.T
	return distance([x0,y0,z0],[x,y,z])

def fitSphere(vertices, p0, nominalRadius, bounds, fitRadius=True):
	if fitRadius:
		errorfun = lambda p,vertices: fittingErrorSphere(p,vertices) - p[3]
	else:
		errorfun = lambda p,vertices: fittingErrorSphere(p,vertices) - nominalRadius

	res = least_squares(errorfun, p0, bounds=([bounds[0][0], bounds[1][0], bounds[2][0], -np.inf], [bounds[0][1], bounds[1][1], bounds[2][1], np.inf]), args=(vertices,))
	centerPoint = res.x[0:3]
	radius = res.x[3]

	return centerPoint, radius

def calculateMeanCurvature(polyData):
	import vtk
	from vtk.util import numpy_support
	cleaner = vtk.vtkCleanPolyData()
	cleaner.SetInputData(polyData)
	cleaner.Update()

	meanCurvature = vtk.vtkCurvatures()
	meanCurvature.SetCurvatureTypeToGaussian()
	meanCurvature.SetInputConnection(cleaner.GetOutputPort())
	meanCurvature.Update()

	npcurv =  numpy_support.vtk_to_numpy(meanCurvature.GetOutput().GetPointData().GetScalars())

	return npcurv

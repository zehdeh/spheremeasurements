import numpy as np
from scipy.optimize import leastsq, least_squares

def getBounds(vertices):
	xbounds = [np.min(vertices[0]),np.max(vertices[0])]
	ybounds = [np.min(vertices[1]),np.max(vertices[1])]
	zbounds = [np.min(vertices[2]),np.max(vertices[2])]
	return [xbounds,ybounds,zbounds]

def distance(p1,p2):
	return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2)

def fittingErrorSphere(center, vertices):
	x0,y0,z0,R = center
	x,y,z = vertices.T
	res = distance([x0,y0,z0],[x,y,z])
	return res

def fitSphere(vertices, nominalRadius, fitRadius=True):
	bounds = getBounds(vertices.T)
	p0 = [bounds[0][0],bounds[1][0],bounds[2][0],nominalRadius]

	if fitRadius:
		errorfun = lambda p,vertices: fittingErrorSphere(p,vertices) - p[3]
	else:
		errorfun = lambda p,vertices: fittingErrorSphere(p,vertices) - nominalRadius

	res = least_squares(errorfun, p0, bounds=([bounds[0][0], bounds[1][0], bounds[2][0], -np.inf], [bounds[0][1], bounds[1][1], bounds[2][1], np.inf]), args=(vertices,))
	centerPoint = res.x[0:3]
	radius = res.x[3]

	return centerPoint, radius

def fittingErrorPlane(point, vertices):
	plane_xyz = point[0:3]
	distance = (plane_xyz*vertices.T).sum(axis=1) + point[3]
	return distance / np.linalg.norm(plane_xyz)

def fitPlane(vertices, p0=[0.1,0.1,0.1,0.1]):
	residuals = lambda p,signal,vertices: fittingErrorPlane(p,vertices)
	res,flag = leastsq(residuals, p0, (None,vertices))
	return res

def calculateCurvature(polyData, ctype=0):
	import vtk
	from vtk.util import numpy_support
	cleaner = vtk.vtkCleanPolyData()
	cleaner.SetInputData(polyData)
	cleaner.Update()

	curvature = vtk.vtkCurvatures()
	if ctype == 0:
		curvature.SetCurvatureTypeToGaussian()
	else:
		curvature.SetCurvatureTypeToMean()
	curvature.SetInputConnection(cleaner.GetOutputPort())
	curvature.Update()

	npcurv =  numpy_support.vtk_to_numpy(curvature.GetOutput().GetPointData().GetScalars())

	return npcurv


def calculateMeanCurvature(polyData):
	return calculateCurvature(polyData, ctype=1)

def calculateGaussianCurvature(polyData):
	return calculateCurvature(polyData, ctype=0)

from shapes import Sphere
from utils import measureRadialDeviation
import numpy as np
from src.fitting import calculateGaussianCurvature, calculateMeanCurvature
from src.spherical_harmonics import getSphericalCoordinates, getCartesianCoordinates
from src.OBJIO import getVTKMesh

class Measure:
	def __init__(self, name, fun, color=False):
		self._name = name
		self.fun = fun
		self._color = color
	@property
	def color(self):
		return self._color
	@property
	def name(self):
		return self._name
	@name.setter
	def name(self, val):
		self._name = val
	def execute(self, sphere):
		res = self.fun(sphere)
		return res

def convertFocusDeviationLabel(label):
	if label[0] == 'p':
		return int(label[1:])
	if label[0] == 'm':
		return int(label[1:]) * -1
	return int(label)

def calculateIdealCurvature(vertices, faces, normals, centerPoint, adaptedRadius):
	phi, theta, r = getSphericalCoordinates(vertices, centerPoint)
	r[:] = adaptedRadius

	vertices = getCartesianCoordinates(phi, theta, r, centerPoint)
	polyData = getVTKMesh(vertices, faces, normals)
	
	return np.std(calculateMeanCurvature(polyData))

def getMeasures(shape):
	measures = []

	fileName = Measure('File', lambda x: x.filePath)
	measures.append(fileName)

	numVertices = Measure('Number of vertices', lambda x: len(x.vertices.T))
	measures.append(numVertices)

	relativeFittingError = Measure('Relative fitting error', 
		lambda x: x.totalFittingError() / x.vertices.shape[1], True)

	#focusDistance = Measure('Focus plane distance', lambda x: int(x.filePath.split('_')[-1][:3]))
	#measures.append(focusDistance)

	#focusDeviation = Measure('Focus plane deviation', lambda x: convertFocusDeviationLabel(x.filePath.split('_')[-3]))
	#measures.append(focusDeviation)

	fittedRadius = Measure('Fitted radius', lambda x: x.fittedRadius)
	measures.append(fittedRadius)

	priorRadius = Measure('Prior radius', lambda x: x.nominalRadius)
	measures.append(priorRadius)

	priorRadius = Measure('Radius difference', lambda x: x.fittedRadius - x.nominalRadius, True)
	measures.append(priorRadius)

	fittingError = Measure('Fitting error total',
		lambda x: x.totalFittingError(), True)
	measures.append(fittingError)

	measures.append(relativeFittingError)

	radialDeviation = Measure('Radial deviation (total)', 
		lambda x: measureRadialDeviation(x.vertices, x.centerPoint, x.fittedRadius)[2], True)
	measures.append(radialDeviation)
	data = []

	idealCurvature = Measure('Ideal Curvature avg',
		lambda x: np.mean(calculateMeanCurvature(x.polyData)) - calculateIdealCurvature(x.vertices.T, x.faces, x.normals, x.centerPoint, x.fittedRadius),True)
	measures.append(idealCurvature)

	avgCurvature = Measure('Mean Curvature avg',
		lambda x: np.mean(calculateMeanCurvature(x.polyData)), True)
	measures.append(avgCurvature)

	stdCurvature = Measure('Mean Curvature std',
		lambda x: np.std(calculateMeanCurvature(x.polyData)), True)
	measures.append(stdCurvature)

	#avgCurvature = Measure('Gaussian Curvature avg',
	#	lambda x: np.mean(calculateGaussianCurvature(x.polyData)), True)
	#measures.append(avgCurvature)

	#stdCurvature = Measure('Gaussian Curvature std',
	#	lambda x: np.std(calculateGaussianCurvature(x.polyData)), True)
	#measures.append(stdCurvature)

	return measures

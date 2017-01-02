from shapes import Sphere
from utils import measureRadialDeviation
import numpy as np
import math
import os
from src.fitting import fittingErrorSphere, calculateGaussianCurvature, calculateMeanCurvature, distance
from src.spherical_harmonics import getSphericalCoordinates, getCartesianCoordinates
from src.OBJIO import getVTKMesh
from src.calibration import getStereoCamerasFromCalibration, StereoCamera, cameraDistance

class Measure:
	def __init__(self, name, fun, color=False, sort=False, formatInteger=False):
		self._name = name
		self.fun = fun
		self._color = color
		self._sort = sort
		self._formatInteger = formatInteger
	@property
	def color(self):
		return self._color
	@property
	def formatInteger(self):
		return self._formatInteger
	@property
	def sort(self):
		return self._sort
	@property
	def name(self):
		return self._name
	@name.setter
	def name(self, val):
		self._name = val
	def execute(self, sphere):
		res = self.fun(sphere)
		return res

def RMSE(x):
	return np.sqrt(np.mean(np.square(distance(x.vertices.T, x.centerPoint) - x.fittedRadius)))

def extractNum(fileName):
	label = fileName[fileName.rfind('_')+1:fileName.find('.')]
	if label[0] == 'P':
		return int(label[1:])*-10
	if label[0] == 'M':
		return int(label[1:]) * 10
	return int(label)*10

def getMeasures(calibration, cameraNo, extract_num):
	measures = []

	if extract_num:
		extracted = Measure('Extracted', lambda x: extractNum(x.fileName), formatInteger=True, sort=True)
		measures.append(extracted)

	numVertices = Measure('Number of vertices', lambda x: len(x.vertices))
	measures.append(numVertices)

	averageRadius = Measure('Average radius',
		lambda x: np.mean(distance(x.vertices.T, x.centerPoint)), False)
	measures.append(averageRadius)

	averageRadiusDifference = Measure('Average radius difference',
		lambda x: x.nominalRadius - np.mean(distance(x.vertices.T, x.centerPoint)), True)
	measures.append(averageRadiusDifference)

	fittedRadius = Measure('Fitted radius', lambda x: x.fittedRadius)
	measures.append(fittedRadius)

	priorRadius = Measure('Prior radius', lambda x: x.nominalRadius)
	measures.append(priorRadius)

	rmse = Measure('RMSE',
		lambda x: RMSE(x), True)
	measures.append(rmse)

	relativeFittingError = Measure('Fitting error mean',
		lambda x: np.mean(np.abs(distance(x.vertices.T, x.centerPoint) - x.fittedRadius)), True)
	measures.append(relativeFittingError)

	relativeFittingErrorStd = Measure('Fitting error std',
		lambda x: np.std(np.fabs(x.fittedRadius - distance(x.vertices.T, x.centerPoint))), True)
	measures.append(relativeFittingErrorStd)

	radialDeviation = Measure('Radial deviation (total)',
		lambda x: np.fabs(np.max(np.fabs(x.nominalRadius - distance(x.vertices.T, x.centerPoint))) - np.min(np.fabs(x.nominalRadius - distance(x.vertices.T, x.centerPoint)))), True)
	measures.append(radialDeviation)
	data = []

	idealCurvature = Measure('Ideal curvature',
		lambda x: 1./x.fittedRadius)
	measures.append(idealCurvature)

	avgCurvature = Measure('Mean curvature avg',
		lambda x: math.fabs(1./x.fittedRadius - np.mean(calculateMeanCurvature(x.polyData))), True)
	measures.append(avgCurvature)

	stdCurvature = Measure('Mean curvature std',
		lambda x: np.std(calculateMeanCurvature(x.polyData)), True)
	measures.append(stdCurvature)

	if calibration is not None and cameraNo is not None:
		cameraFocusCalibrationPath = calibration
		stereoCameras = getStereoCamerasFromCalibration(cameraFocusCalibrationPath)
		cameraPosition = (stereoCameras[cameraNo].A.position + stereoCameras[cameraNo].B.position) / 2

		baseLine = distance(stereoCameras[cameraNo].A.position, stereoCameras[cameraNo].B.position)

		measures.append(Measure('Distance from camera', lambda x: distance(x.centerPoint,cameraPosition), False, True))

	return measures

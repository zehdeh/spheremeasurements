from shapes import Sphere, Plane
from utils import measureRadialDeviation

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

def getMeasures(shape):
	measures = []

	fileName = Measure('File', lambda x: x.filePath)
	measures.append(fileName)

	numVertices = Measure('Number of vertices', lambda x: len(x.vertices.T))
	measures.append(numVertices)

	relativeFittingError = Measure('Relative fitting error', 
		lambda x: x.totalFittingError() / x.vertices.shape[1], True)

	if shape is Sphere:
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

		fittingError = Measure('Fitting error', 
			lambda x: x.totalFittingError())
		measures.append(fittingError)

		measures.append(relativeFittingError)

		radialDeviation = Measure('Radial deviation (total)', 
			lambda x: measureRadialDeviation(x.vertices, x.centerPoint, x.fittedRadius)[2], True)
		measures.append(radialDeviation)
		data = []

		curvature = Measure(['Average curvature', 'Standard deviation'], 
			lambda x: x.measureCurvature(), True)
		measures.append(curvature)
	if shape is Plane:
		fittingError = Measure('Fitting error', 
			lambda x: x.totalFittingError(), True)
		measures.append(fittingError)

		fittingError = Measure('Plane std deviation', 
			lambda x: x.planeDeviation(), True)
		measures.append(fittingError)

		measures.append(relativeFittingError)

		curvature = Measure(['Average curvature', 'Standard deviation'], 
			lambda x: x.measureCurvature(), True)
		measures.append(curvature)

	return measures

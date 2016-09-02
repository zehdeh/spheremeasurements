import numpy as np
from scipy.optimize import leastsq, least_squares
from opendr.serialization import load_mesh
from utils import distance
from opendr.geometry import GaussianCurvature

def fitfunc(center, coords):
	x0,y0,z0,R = center
	x,y,z = coords
	return distance([x0,y0,z0],[x,y,z]) - R
def radiusErrorFunction(p,x): return fitfunc(p,x)

class Shape(object):
	def __init__(self,filePath):
		self._mesh = load_mesh(filePath)
		self._vertices = self.mesh.v.T
		self._faces = self.mesh.f
		self._filePath = filePath
	@property
	def mesh(self):
		return self._mesh
	@mesh.setter
	def mesh(self,value):
		self._mesh = value
	@property
	def vertices(self):
		return self._vertices
	@mesh.setter
	def vertices(self,value):
		self._vertices = value
	@property
	def filePath(self):
		return self._filePath
	@filePath.setter
	def filePath(self,value):
		self_filePath = value

def residuals(p,signal,vertices):
	return fitPlaneFun(vertices,p)

def fitPlaneFun(vertices,p):
	plane_xyz = p[0:3]
	distance = (plane_xyz*vertices.T).sum(axis=1) + p[3]
	return distance / np.linalg.norm(plane_xyz)

class Plane(Shape):
	def __init__(self,filePath):
		super(Plane,self).__init__(filePath)

		res = self.fitPlane()
		print res
	@property
	def vertices(self):
		return self._vertices
	def fitPlane(self):
		p0 = [0.1,0.1,0.1,0.1]
		res = leastsq(residuals, p0, (None,self._vertices))[0]
		return res

class Sphere(Shape):
	def __init__(self, filePath, nominalRadius):
		Shape.__init__(self,filePath)

		fittedRadius, centerPoint = self.fitSphere()
		self._fittedRadius = fittedRadius
		self._centerPoint = centerPoint
		self._nominalRadius = nominalRadius
		self.curvature = np.asarray(GaussianCurvature(self._vertices.T, self._faces))
	@property
	def vertices(self):
		return self._vertices
	@property
	def fittedRadius(self):
		return self._fittedRadius
	@fittedRadius.setter
	def fittedRadius(self,value):
		self._fittedRadius = value
	@property
	def centerPoint(self):
		return self._centerPoint
	@centerPoint.setter
	def centerPoint(self,value):
		self._centerPoint = value
	@property
	def nominalRadius(self):
		return self._nominalRadius
	@nominalRadius.setter
	def nominalRadius(self, value):
		self._nominalRadius = value
	def getPointBounds(self):
		xbounds = [np.min(self._vertices[0]),np.max(self._vertices[0])]
		ybounds = [np.min(self._vertices[1]),np.max(self._vertices[1])]
		zbounds = [np.min(self._vertices[2]),np.max(self._vertices[2])]
		return [xbounds,ybounds,zbounds]
	def measureCurvature(self):
		return [
			np.mean(self.curvature),
			np.std(self.curvature)]
	def totalFittingError(self):
		return np.sum(np.abs(self._fittedRadius - distance(self._vertices,self._centerPoint)))
		#return np.sum(np.abs(radiusErrorFunction(self.centerPoint.tolist() + [self.fittedRadius], self._vertices)))/len(self._vertices.T)
	def fitSphere(self):
		#vertices = generateSphere(300, 35)
		pointBounds = self.getPointBounds()

		centerPointGuess = [pointBounds[0][0],pointBounds[1][0],pointBounds[2][0],36]
		#centerPointGuess = [-124.63678821,-283.53124005,-334.92304383, 1]
		#res = least_squares(errfunc, centerPointGuess, bounds=([pointBounds[0][0],pointBounds[1][0],pointBounds[2][0],1],
		#[pointBounds[0][1],pointBounds[1][1],pointBounds[2][1],np.inf]), args=(vertices,))
		res = least_squares(radiusErrorFunction, centerPointGuess, bounds=(-np.inf, np.inf), args=(self._vertices,))
		centerPoint = res.x

		fittedRadius = centerPoint[3]
		centerPoint = centerPoint[0:3]

		return fittedRadius,centerPoint

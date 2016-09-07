import abc
import numpy as np
from scipy.optimize import leastsq, least_squares
from opendr.serialization import load_mesh
from utils import distance
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from opendr.geometry import GaussianCurvature
from mayavi.mlab import *

def randomPartition(n, nData):
	allIdxs = np.arange(nData)
	np.random.shuffle(allIdxs)
	idxs1 = allIdxs[:n]
	idxs2 = allIdxs[n:]
	return idxs1, idxs2

def ransac(data,fitfun,errorfun,n,k,t,d):
	debug = False
	iterations = 0
	bestfit = None
	besterr = np.inf
	best_inlier_idxs= None
	while iterations < k:
		maybe_idxs, test_idxs = randomPartition(n,data.shape[0])
		maybeinliers = data[maybe_idxs,:]
		test_points = data[test_idxs]
		maybemodel = fitfun(maybeinliers)
		test_err = errorfun(maybemodel, test_points)
		also_idxs = test_idxs[test_err < t] # select indices of rows with accepted points
		alsoinliers = data[also_idxs,:]
		if debug:
			print 'test_err.min()',test_err.min()
			print 'test_err.max()',test_err.max()
			print 'numpy.mean(test_err)',np.mean(test_err)
			print 'iteration %d:len(alsoinliers) = %d'%(
				iterations,len(alsoinliers))
		if len(alsoinliers) > d:
			betterdata = np.concatenate( (maybeinliers, alsoinliers) )
			bettermodel = fitfun(betterdata)
			better_errs = errorfun(bettermodel, betterdata)
			thiserr = np.mean( better_errs )
			if thiserr < besterr:
				bestfit = bettermodel
				besterr = thiserr
				best_inlier_idxs = np.concatenate( (maybe_idxs, also_idxs) )
		iterations+=1
	if bestfit is None:
		raise ValueError("did not meet fit acceptance criteria")
	else:
		return bestfit

class Shape(object):
	__metaclass__ = abc.ABCMeta
	def __init__(self,filePath):
		self._mesh = load_mesh(filePath)
		self._vertices = self.mesh.v.T
		self._faces = self.mesh.f
		self._filePath = filePath

		self.curvature = np.asarray(GaussianCurvature(self._vertices.T, self._faces))
	@property
	def mesh(self):
		return self._mesh
	@mesh.setter
	def mesh(self,value):
		self._mesh = value
	@property
	def faces(self):
		return self._faces
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
	@abc.abstractmethod
	def totalFittingError(self):
		return
	def render(self,fig):
		from mayavi import mlab

		mesh = triangular_mesh(self._vertices[0], self._vertices[1], self._vertices[2], self._faces, scalars=self.curvature)
		mlab.scalarbar(mesh)
		mlab.outline()

		mlab.show()
		"""
		means = [np.mean(v) for v in self.vertices]
		relativeMins = [np.min(v - m) for v,m in zip(self.vertices,means)]
		relativeMaxs = [np.max(v - m) for v,m in zip(self.vertices,means)]
		limitsMin = np.add(means,relativeMins)
		limitsMax = np.add(means,relativeMaxs)
		valMin = np.min(np.add(means,relativeMins))
		valMax = np.max(np.add(means,relativeMaxs))

		poly3dCollection = Poly3DCollection(self._vertices.T[self._faces], edgecolors='black')
		ax3d = fig.add_subplot(111, projection='3d')
		
		ax3d.add_collection(poly3dCollection)

		ax3d.set_xlim(limitsMin[0], limitsMax[0])
		ax3d.set_ylim(limitsMin[1], limitsMax[1])
		ax3d.set_zlim(limitsMin[2], limitsMax[2])
		

		ax3d.set_xlabel('x')
		ax3d.set_ylabel('y')
		ax3d.set_zlabel('z')

		return ax3d
		"""

class Plane(Shape):
	def __init__(self,filePath):
		super(Plane,self).__init__(filePath)

		self.model = np.asarray(self.fitPlane())
	@property
	def vertices(self):
		return self._vertices
	def pointDistance(self, point, vertices):
		plane_xyz = point[0:3]
		distance = (plane_xyz*vertices.T).sum(axis=1) + point[3]
		return distance / np.linalg.norm(plane_xyz)
	def fitPlane(self, doRansac = False):
		p0 = [0.1,0.1,0.1,0.1]
		errorfun = lambda p,vertices: self.pointDistance(p,vertices)
		residuals = lambda p,signals,vertices: self.pointDistance(p,vertices)
		fitfun = lambda data: leastsq(residuals, p0, (None,self._vertices))[0]
		if doRansac:
			res = ransac(self._vertices.T, fitfun, errorfun, int(0.1*len(self._vertices.T)), 100, 0.01, int(0.5*len(self._vertices.T)))
		else:
			res = fitfun(self._vertices)
		return res
	def totalFittingError(self):
		return self.pointDistance(self.model, self._vertices).sum()
	def planeDeviation(self):
		maxDeviation = np.max(self.pointDistance(self.model, self._vertices))
		minDeviation = np.min(self.pointDistance(self.model, self._vertices))
		return abs(maxDeviation-minDeviation)
	def measureCurvature(self):
		return [
			np.mean(self.curvature),
			np.std(self.curvature)]
	def render(self,fig):
		super(Plane,self).render(fig)
		"""
		normal = self.model[:3]
		d = -np.array([0,0,0]).dot(normal)#d = self.model[3]
		x = np.linspace(-1,1,5)
		y = np.linspace(-1,1,5)
		xx,yy = np.meshgrid(x,y)
		z = ((-normal[0]*xx - normal[1]*yy - d)*1./normal[2])
		yy = yy

		ax3d.plot_wireframe(xx,yy,z)
		"""

class Sphere(Shape):
	def __init__(self, filePath, nominalRadius):
		Shape.__init__(self,filePath)

		fittedRadius, centerPoint = self.fitSphere(True)
		self._fittedRadius = fittedRadius
		self._centerPoint = centerPoint
		self._nominalRadius = nominalRadius
		print self.totalFittingError()
		fittedRadius, centerPoint = self.fitSphere(False)
		self._fittedRadius = fittedRadius
		self._centerPoint = centerPoint
		self._nominalRadius = nominalRadius
		print self.totalFittingError()

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
	def fittingError(self,center, vertices):
		x0,y0,z0,R = center
		x,y,z = vertices.T
		return distance([x0,y0,z0],[x,y,z]) - R
	def totalFittingError(self):
		return np.sum(np.abs(self._fittedRadius - distance(self._vertices,self._centerPoint)))
		#return np.sum(np.abs(radiusErrorFunction(self.centerPoint.tolist() + [self.fittedRadius], self._vertices)))/len(self._vertices.T)
	def fitSphere(self, doRansac = False):
		#vertices = generateSphere(300, 35)
		pointBounds = self.getPointBounds()
		errorfun = lambda p,vertices: self.fittingError(p,vertices)

		centerPointGuess = [pointBounds[0][0],pointBounds[1][0],pointBounds[2][0],36]
		#centerPointGuess = [-124.63678821,-283.53124005,-334.92304383, 1]
		#res = least_squares(errfunc, centerPointGuess, bounds=([pointBounds[0][0],pointBounds[1][0],pointBounds[2][0],1],
		#[pointBounds[0][1],pointBounds[1][1],pointBounds[2][1],np.inf]), args=(vertices,))
		fitfun = lambda data: least_squares(errorfun, centerPointGuess, bounds=(-np.inf, np.inf), args=(data,)).x
		if doRansac:
			res = ransac(self._vertices.T, fitfun, errorfun, int(0.1*len(self._vertices.T)), 100, 0.01, int(0.5*len(self._vertices.T)))
		else:
			res = fitfun(self._vertices.T)

		centerPoint = res

		fittedRadius = centerPoint[3]
		centerPoint = centerPoint[0:3]

		return fittedRadius,centerPoint

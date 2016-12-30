import abc
import numpy as np
from scipy.optimize import leastsq, least_squares
from opendr.serialization import load_mesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from opendr.geometry import GaussianCurvature
import vtk
from src.fitting import distance, calculateMeanCurvature, fitSphere
from src.OBJIO import loadOBJviaVTK, loadOBJ, getVTKMesh


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

class Sphere(object):
	def __init__(self, filePath, nominalRadius, fitRadius = True):
		self._filePath = filePath

		self._vertices, self._faces, self._normals, self.polyData = loadOBJviaVTK(filePath)
		#self._vertices, self._faces, self._normals = loadOBJ(filePath)
		#self.polyData = getVTKMesh(self._vertices, self._faces, self._normals)

		self.curvature = calculateMeanCurvature(self.polyData)

		self._nominalRadius = float(nominalRadius)
		centerPoint, fittedRadius = fitSphere(self._vertices, nominalRadius, False)
		self._fittedRadius = fittedRadius
		self._centerPoint = centerPoint
	@property
	def faces(self):
		return self._faces
	@property
	def normals(self):
		return self._normals
	@property
	def filePath(self):
		return self._filePath
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
		return distance([x0,y0,z0],[x,y,z])
	def totalFittingError(self):
		return np.sum(np.fabs(self._fittedRadius - distance(self._vertices,self._centerPoint)))
		#return np.sum(np.abs(radiusErrorFunction(self.centerPoint.tolist() + [self.fittedRadius], self._vertices)))/len(self._vertices.T)

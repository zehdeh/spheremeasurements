from math import sin,cos,sqrt,acos,fabs
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from src.fitting import distance
import numpy as np
import scipy

def rodrigues(src):
	dst = np.zeros(3)
	matJ = np.zeros((3,9))

	matR = np.zeros((3,3))
	matU = np.zeros((3,3))
	matV = np.zeros((3,3))
	matW = np.zeros(3)

	matR = src
	matU,matW,matV = np.linalg.svd(matR)
	#print 'U:'
	#print matU
	#print 'V:'
	#print matV
	matR = matU.dot(matV.T)

	R = matR.ravel()
	rx = R[7] - R[5]
	ry = R[2] - R[6]
	rz = R[3] - R[1]

	s = sqrt((rx*rx + ry*ry + rz*rz)*0.25)
	c = (R[0] + R[4] + R[8] - 1)*0.5
	c = 1. if c > 1. else (-1. if c < -1. else c)
	theta = acos(c)

	if s < 1e-5:
		if c > 0:
			rx = 0
			ry = 0
			rz = 0
		else:
			t = (R[0] + 1)*0.5
			rx = sqrt(max(t,0))
			t = (R[4] + 1)*0.5
			ry = sqrt(max(t,0))*(-1. if R[1] < 0 else 1.)
			t = (R[8] + 1)*0.5
			rz = sqrt(max(t,0))*(-1. if R[2] < 0 else 1.)
			if fabs(rx) < fabs(ry) and fabs(rx) < fabs(rz) and (R[5] > 0) != (ry*rz > 0):
				rz = -rz
			theta /= sqrt(rx*rx + ry*ry + rz*rz)
			rx *= theta
			ry *= theta
			rz *= theta
	else:
		vth = 1./(2*s)
		vth *= theta

		rx *= vth
		ry *= vth
		rz *= vth
	
	dst[0] = rx
	dst[1] = ry
	dst[2] = rz

	return dst


def cartesianProduct(arrays, out=None):
	arrays = [np.asarray(x) for x in arrays]
	dtype = arrays[0].dtype

	n = np.prod([x.size for x in arrays])
	if out is None:
		out = np.zeros([n, len(arrays)], dtype=dtype)

	m = n / arrays[0].size
	out[:,0] = np.repeat(arrays[0], m)
	if arrays[1:]:
		cartesianProduct(arrays[1:], out=out[0:m,1:])
		for j in xrange(1, arrays[0].size):
			out[j*m:(j+1)*m,1:] = out[0:m,1:]
	return out.astype(int)

def scanCircleFit(centerPoints):
	def fitFuncSphere(p, coords):
		x0,y0,z0,R = p
		x,y,z = coords.T
		return distance([x0,y0,z0],[x,y,z]) - R
	
	errorFunc = lambda p,coords: fitFuncSphere(p,coords)

	p0 = [0,0,0,2000]
	bounds = ([-np.inf,-np.inf,-np.inf,2000],[np.inf,np.inf,np.inf,np.inf])
	res = scipy.optimize.least_squares(errorFunc, p0, bounds=bounds,args=(centerPoints,))
	print res.x
	return res.x

def scanLineFit(centerPoints):
	centerPoints = np.asarray(centerPoints).T
	mean = centerPoints.mean(axis=1)
	uu,dd,vv = np.linalg.svd(centerPoints.T - mean)
	AP = centerPoints.T - mean
	pointDistances = AP.dot(vv[0])[..., None]
	projectedPoints = ((pointDistances * vv[0][None,...]) + mean).T

	return mean, vv[0], projectedPoints, pointDistances

def toSphericalCoordinates(xyz):
	ptsnew = np.zeros(xyz.shape)
	xy = xyz[:,0]**2 + xyz[:,1]**2
	ptsnew[:,0] = np.sqrt(xy + xyz[:,2]**2)
	ptsnew[:,1] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
	#ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
	ptsnew[:,2] = np.arctan2(xyz[:,1], xyz[:,0])
	return ptsnew


def measureRadialDeviation(sphere, center, radius):
	minRadius = np.inf
	maxRadius = -np.inf

	for i in sphere.T:
		currentRadius = distance(i,center)

		if minRadius > currentRadius:
			minRadius = currentRadius

		if maxRadius < currentRadius:
			maxRadius = currentRadius
	
	return [abs(radius - minRadius), abs(radius - maxRadius), abs(maxRadius - minRadius)]


class Arrow3D(FancyArrowPatch):
	def __init__(self, xs, ys, zs, *args, **kwargs):
		FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
		self._verts3d = xs, ys, zs

	def draw(self, renderer):
		xs3d, ys3d, zs3d = self._verts3d
		xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
		self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
		FancyArrowPatch.draw(self, renderer)

def generateSphere(samples=1,radius=1,randomize=True):
	rnd = 1.
	if randomize:
		rnd = random.random() * samples
	
	points = []
	offset = 2./samples
	increment = math.pi * (3. - math.sqrt(5.))

	for i in range(samples):
#		if(((i * offset) - 1) + (offset / 2)) < 0:
#			continue
		y = ((i * offset) - 1) + (offset / 2)
		r = math.sqrt(1 - pow(y,2))

		phi = ((i + rnd) % samples) * increment

		y *= radius
		x = math.cos(phi) * r * radius
		z = math.sin(phi) * r * radius
		points.append([x,y,z])
	vertices = np.array(points)

	cvx = ConvexHull(vertices)
	faces = cvx.simplices
	return vertices.T, faces



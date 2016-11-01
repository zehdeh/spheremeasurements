from math import sin,cos
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from src.fitting import distance
import numpy as np
import scipy

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



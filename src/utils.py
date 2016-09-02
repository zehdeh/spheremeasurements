from matplotlib.patches import FancyArrowPatch
import numpy as np
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

def distance(p1,p2):
	return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2)


#! /usr/bin/python

#from mayavi.mlab import *
import sys
import math, random
from mpl_toolkits.mplot3d import Axes3D, proj3d
import matplotlib.pyplot as plt
from matplotlib.mlab import PCA
from matplotlib.patches import FancyArrowPatch
from matplotlib.tri import Triangulation
import numpy as np
import numpy.linalg
from scipy.optimize import leastsq, least_squares
from scipy.spatial import Delaunay, ConvexHull
import pymesh

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
	
	return np.array(points).T

def distance(p1,p2):
	return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2)

def loadOBJ(filename):
	numVerts = 0
	verts = []
	norms = []
	vertsOut = []
	normsOut = []
	for line in open(filename, "r"):
		vals = line.split()
		if len(vals) == 0: continue
		if vals[0] == "v":
			v = map(float, vals[1:4])
			verts.append(v)
		if vals[0] == "vn":
			n = map(float, vals[1:4])
			norms.append(n)
		if vals[0] == "f":
			for f in vals[1:]:
				w = f.split("/")
				# OBJ Files are 1-indexed so we must subtract 1 below
				vertsOut.append(list(verts[int(w[0])-1]))
				normsOut.append(list(norms[int(w[2])-1]))
				numVerts += 1
	return np.array(verts).T, np.array(norms).T

def getPointBounds(coords):
	xbounds = [np.min(coords[0]),np.max(coords[0])]
	ybounds = [np.min(coords[1]),np.max(coords[1])]
	zbounds = [np.min(coords[2]),np.max(coords[2])]
	return [xbounds,ybounds,zbounds]

def measureDiameter(sphere):
	#vertices = generateSphere(300, 35)
	pointBounds = getPointBounds(sphere)

	centerPointGuess = [pointBounds[0][0],pointBounds[1][0],pointBounds[2][0],36]
	#centerPointGuess = [-124.63678821,-283.53124005,-334.92304383, 1]
	errfunc = lambda p,x: fitfunc(p,x)
	#res = least_squares(errfunc, centerPointGuess, bounds=([pointBounds[0][0],pointBounds[1][0],pointBounds[2][0],1],
	#[pointBounds[0][1],pointBounds[1][1],pointBounds[2][1],np.inf]), args=(vertices,))
	res = least_squares(errfunc, centerPointGuess, bounds=(-np.inf, np.inf), args=(sphere,))
	centerPoint = res.x

	lsqerr = np.sum(errfunc(centerPoint,sphere)**2)
	noPoints = np.size(sphere)

	print "Nominal error (least squares): " + str(lsqerr)
	print "Number of points: " + str(noPoints)
	print "Relative error: " + str(lsqerr/noPoints)
	print "D_A: " + str(centerPoint[3]*2)
	return centerPoint

def measureRadialDeviation(sphere, center, radius):
	minRadius = np.inf
	maxRadius = -np.inf

	for i in sphere.T:
		currentRadius = distance(i,center)

		if minRadius > currentRadius:
			minRadius = currentRadius

		if maxRadius < currentRadius:
			maxRadius = currentRadius
	
	return [abs(radius - minRadius), abs(radius - maxRadius)]

def plotRadialDeviation(fig,sphere, center, radius, minRadius, maxRadius):
	zValues = sphere[2]-center[2]
	h = np.sqrt((sphere[0]-center[0])**2 + (sphere[2]-center[2])**2)
	r = distance(sphere,center)

	alphas = np.arcsin(zValues/h)
	zIndices = np.where(zValues<0)
	z2Indices = np.where(zValues>0)
	
	xIndices = np.where((sphere[0]-center[0])<0)
	#alphas[zIndices] += np.pi
	alphas[xIndices] = alphas[xIndices]*-1+np.pi

	ax = fig.add_subplot(211, projection='polar')
	ax.scatter(alphas, r)
	ax.grid(True)

def fitfunc(center, coords):
	x0,y0,z0,R = center
	x,y,z = coords
	return distance([x0,y0,z0],[x,y,z]) - R

if __name__ == '__main__':
	if len(sys.argv) < 2:
		#print "Please specify a file to open"
		sys.exit(1)

	#mesh = pymesh.load_mesh(sys.argv[1])
	#vertices = mesh.vertices.T
	#print mesh.faces

	#vertices,norms = loadOBJ(sys.argv[1])
	vertices = generateSphere(300, 35)

	cvx = ConvexHull(vertices.T)
	x,y,z = vertices

	tri = Triangulation(x,y, triangles=cvx.simplices)
	#tri = Delaunay(vertices.T)
	mesh = pymesh.form_mesh(vertices, tri.triangles)
	vertices = mesh.vertices
	#print tri.simplices[1,:]

	mesh.add_attribute("vertex_gaussian_curvature")
	#curvature = mesh.get_attribute("vertex_gaussian_curvature")
	#print curvature

	cov = np.cov(vertices)
	eigval, eigvec = np.linalg.eig(cov)

	centerPoint = measureDiameter(vertices)
	meanRadius = centerPoint[3]
	[minRadius,maxRadius] = measureRadialDeviation(vertices, centerPoint, meanRadius)

	print 'Expected radius: ' + str(35)
	print 'Fitted radius: ' + str(meanRadius)
	print 'Radial deviation: (+ ' + str(maxRadius) + ',- ' + str(minRadius) + ')'

	fig = plt.figure()
	plotRadialDeviation(fig, vertices, centerPoint, meanRadius, minRadius, maxRadius)
	plt.axis('scaled')
	ax = fig.add_subplot(212, projection='3d')
	#ax.scatter(vertices[0].tolist(), vertices[1].tolist(), vertices[2].tolist(), c="b")
	ax.plot_trisurf(vertices[0].tolist(), vertices[1].tolist(), vertices[2].tolist(), triangles=tri.triangles)
	ax.scatter(centerPoint[0], centerPoint[1], centerPoint[2], c='black')

	for vec,val in zip(eigvec.T,eigval):
		arrowEndPoint = centerPoint[0:3] + vec*val*0.1
		arrow = Arrow3D([centerPoint[0], arrowEndPoint[0]], [centerPoint[1], arrowEndPoint[1]], [centerPoint[2], arrowEndPoint[2]], mutation_scale=20, lw=2, arrowstyle='-|>', color="r")
		ax.add_artist(arrow)

	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')
	plt.show()

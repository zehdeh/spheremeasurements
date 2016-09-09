import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from src.mplcolorhelper import MplColorHelper
from src.shapes import Sphere,Plane

def getCurvatureColor(shape, normalize=False):
	maxCurvature = np.mean(shape.curvature[shape.faces], axis=1)
	if normalize:
		maxCurvature = maxCurvature / np.max(maxCurvature)
	colorHelper = MplColorHelper('coolwarm', np.min(maxCurvature), np.max(maxCurvature), maxCurvature)
	polyColors = [colorHelper.get_rgb(c) for c in maxCurvature]
	return polyColors

def plotPolygons(shape):
	#fig = plt.figure()
	print np.mean(shape.curvature)
	print np.std(shape.curvature)
	shape.render(1)

	"""
	curv = shape.curvature
	a = int(math.sqrt(shape.vertices.shape[1]))
	print a
	curv = np.reshape(curv, (a,a))

	S = np.fft.rfft(curv)
	F2 = np.fft.fftshift(S)

	psf2D = np.abs(F2)**2

	fig, ax = plt.subplots()
	plt.imshow(np.log10(psf2D))
	plt.show()
	#plt.show()
	"""

if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please provide as an argument the file you want to plot'
		sys.exit(1)
	shapeObj = Plane(sys.argv[1])
	plotPolygons(shapeObj)

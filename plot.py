#! /usr/bin/python

import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from src.mplcolorhelper import MplColorHelper
from src.shapes import Sphere,Plane

def plotPolygons(shape):
	#fig = plt.figure()
	shape.render()

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
	shapeObj = Sphere(sys.argv[1],300)
	plotPolygons(shapeObj)

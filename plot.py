#! /usr/bin/python

import sys
import numpy as np
import math
from src.fitting import distance
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from src.mplcolorhelper import MplColorHelper
from src.shapes import Sphere,Plane


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please provide as an argument the file you want to plot'
		sys.exit(1)
	shape = Sphere(sys.argv[1],150)
	if len(sys.argv) == 2:
		sys.argv.append('error')

	if sys.argv[2] == 'curvature':
		curvature = shape.getCurvature()
		plt.hist(curvature,bins=50,normed=1)
		plt.title("Curvature")
		plt.show()
		shape.render(curvature)
	elif sys.argv[2] == 'gradient':
		curvature = shape.getCurvature()
		grad = np.gradient(curvature)
		shape.render(grad)
	else:
		error = np.abs(shape.fittedRadius - distance(shape.vertices,shape.centerPoint))
		plt.hist(error,bins=50,normed=1)
		plt.title("Error")
		plt.show()
		shape.render(error)

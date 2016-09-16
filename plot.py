#! /usr/bin/python

import sys
import numpy as np
import math
from src.utils import distance
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from src.mplcolorhelper import MplColorHelper
from src.shapes import Sphere,Plane


if __name__ == '__main__':
	if len(sys.argv) < 2:
		print 'Please provide as an argument the file you want to plot'
		sys.exit(1)
	shape = Sphere(sys.argv[1],150)
	if sys.argv[2] == 'curvature':
		shape.render(shape.getCurvature())
	else:
		error = np.abs(shape.fittedRadius - distance(shape.vertices,shape.centerPoint))
		shape.render(error)

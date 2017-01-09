#! /usr/bin/python

import numpy as np
import argparse
from matplotlib import pyplot as plt

from src.utils import checkFile

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='plots a 3d grid in 2d')
	parser.add_argument("file", help="The npy file to plot", type=lambda x: checkFile(x,'.npy'))
	parser.add_argument("axis", help="The axis in which you want to reduce the grid", type=int)
	parser.add_argument("--verbose", help="Show debug information", action='store_true')
	args = parser.parse_args()

	if args.axis < 0 or args.axis > 3:
		raise argparse.ArgumentTypeError("Invalid axis: {0}".format(args.axis))

	dataMatrix = np.load(args.file)
	print np.unravel_index(np.argmax(dataMatrix),[50,50,50])
	print np.max(dataMatrix)
	print np.mean(dataMatrix)
	twoDimMatrix = np.sum(dataMatrix,axis=args.axis).astype(np.float)
	twoDimMatrix /= twoDimMatrix.max()

	fig = plt.figure(facecolor='white')
	if args.axis == 0:
		temp = twoDimMatrix[::-1,:]
		twoDimMatrix[:,:] = temp
		plt.xlabel('z')
		plt.ylabel('y')
	elif args.axis == 1:
		plt.xlabel('x')
		plt.ylabel('z')
	else:
		twoDimMatrix = twoDimMatrix.T
		temp = twoDimMatrix[::-1,:]
		twoDimMatrix[:,:] = temp
		plt.xlabel('x')
		plt.ylabel('y')
	plt.imshow(twoDimMatrix, interpolation='none')
	plt.colorbar()
	plt.show()


import numpy as np

class Mesh(object):
	def __init__(self, vertices, faces, normals):
		self.vertices = vertices
		self.faces = faces
		self.normals = normals
	def getBounds(self):
		xbounds = [np.min(self.vertices[0]),np.max(self.vertices[0])]
		ybounds = [np.min(self.vertices[1]),np.max(self.vertices[1])]
		zbounds = [np.min(self.vertices[2]),np.max(self.vertices[2])]
		return [xbounds,ybounds,zbounds]

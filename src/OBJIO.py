from array import *
import vtk
import locale
from vtk.util import numpy_support
import numpy as np

def loadOBJviaVTK(fileName):

	locale.setlocale(locale.LC_NUMERIC, 'C')
	reader = vtk.vtkOBJReader()
	reader.SetFileName(fileName)
	reader.Update()

	polyData = reader.GetOutput()

	vertices =  numpy_support.vtk_to_numpy(polyData.GetPoints().GetData())
	faces = numpy_support.vtk_to_numpy(polyData.GetPolys().GetData())
	faces = faces.reshape((-1,4)).T[1:4].T
	normals = numpy_support.vtk_to_numpy(polyData.GetPointData().GetNormals())

	return vertices, faces, normals, polyData

def loadOBJ(fileName):
	"""Loads a Wavefront OBJ file. """
	vertices = []
	normals = []
	texcoords = []
	faces = []

	material = None
	for line in open(fileName, "r"):
		if line.startswith('#'): continue
		values = line.split()
		if not values: continue
		if values[0] == 'v':
			v = list(map(float, values[1:4]))
			vertices.append(v)
		elif values[0] == 'vn':
			v = list(map(float, values[1:4]))
			normals.append(v)
		elif values[0] == 'vt':
			texcoords.append((map(float, values[1:3])))
		elif values[0] == 'f':
			face = []
			texcoords = []
			norms = []
			for v in values[1:]:
				if '//' in v:
					glue = '//'
				else:
					glue = '/'
				w = v.split(glue)
				face.append(int(w[0]) - 1)
				if len(w) >= 2 and len(w[1]) > 0:
					texcoords.append(int(w[1]))
				else:
					texcoords.append(0)
					if len(w) >= 3 and len(w[2]) > 0:
						norms.append(int(w[2]))
					else:
						norms.append(0)
			faces.append(face)
	return np.asarray(vertices), np.asarray(faces), np.asarray(normals)

def writeOBJ(fileName, vertices, faces, normals):
	with open(fileName, 'w') as f:
		f.write("# OBJ file\n")
		for v in vertices:
			f.write('v ' +' '.join([format(x,'.4f') for x in v]))
			f.write("\n")
		for vn in normals:
			f.write('vn ' +' '.join([format(x,'.4f') for x in vn]))
			f.write("\n")
		for p in faces:
			f.write("f")
			for i in p:
				f.write(" %d" % (i + 1))
			f.write("\n")

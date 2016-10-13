#!/usr/bin/env python
 
import sys
import vtk
import numpy as np
from vtk.util import numpy_support
from src.fitting import fitSphere, fittingErrorSphere
from src.mesh import Mesh
from src.OBJIO import loadOBJ, loadOBJviaVTK
import matplotlib.pyplot as plt
from pymetis import part_graph
from opendr.topology import get_vert_connectivity
 
class CurvaturesDemo():
	def CurvaturesDemo(self):
		#vertices, faces, normals = loadOBJ(sys.argv[1])
		vertices, faces, normals, polyData = loadOBJviaVTK(sys.argv[1])

		bounds = Mesh(vertices.T, faces, normals).getBounds()
		p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
		cp, radius = fitSphere(vertices,p0,150, bounds)
		errors = fittingErrorSphere(cp.tolist() + [radius], vertices) - 150

		reader = vtk.vtkOBJReader()
		reader.SetFileName(sys.argv[1])


		reader.Update()
		polydata = reader.GetOutput()

		# Colour transfer function.
		ctf = vtk.vtkColorTransferFunction()
		ctf.SetColorSpaceToDiverging()
		ctf.AddRGBPoint(0.0, 0.230, 0.299, 0.754)
		ctf.AddRGBPoint(0.5, 1.0, 1.0, 1.0)
		ctf.AddRGBPoint(1.0, 0.706, 0.016, 0.150)
		cc = list()
		for i in range(256):
			cc.append(ctf.GetColor(float(i) / 255.0)) 

		scalars = vtk.vtkFloatArray()
		scalars.SetNumberOfComponents(1)

		for errorVal in errors:
			scalars.InsertNextTupleValue([errorVal])

		mean = np.mean(errors)
		std = np.std(errors)
		min3 = mean - std
		max3 = mean + std

		polydata2 = vtk.vtkPolyData()
		polydata2.DeepCopy(polydata)

		cnct = get_vert_connectivity(vertices, faces)
		nbrs = [np.nonzero(np.array(cnct[:,i].todense()))[0] for i in range(cnct.shape[1])]
		edgeCuts, parts = part_graph(5,nbrs)

		partsvtk = vtk.vtkUnsignedCharArray()
		partsvtk.SetNumberOfComponents(1)
		for p in parts:
			partsvtk.InsertNextTuple([p])

		polydata.GetPointData().SetScalars(scalars)
		polydata2.GetPointData().SetScalars(partsvtk)

		# Now we have the sources, lets put them into a list.
		sources = list()
		sources.append(reader)
		sources.append(reader)
		sources.append(polydata)
		sources.append(polydata2)

		curvatures = list()        
		for idx in range(2):
			curvatures.append(vtk.vtkCurvatures())
			curvatures[idx].SetInputConnection(sources[idx].GetOutputPort())
			if idx % 2 == 0:
				curvatures[idx].SetCurvatureTypeToGaussian()
				curvatures[idx].Update()
				npcurv1 =  numpy_support.vtk_to_numpy(curvatures[idx].GetOutput().GetPointData().GetScalars())
		
				mean = np.mean(npcurv1)
				std = np.std(npcurv1)
				print std

				min1 = 0.01*(mean - std)
				max1 = 0.01*(mean + std)
			else:
				curvatures[idx].SetCurvatureTypeToMean()
				curvatures[idx].Update()
				npcurv2 =  numpy_support.vtk_to_numpy(curvatures[idx].GetOutput().GetPointData().GetScalars())


				mean = np.mean(npcurv2)
				std = np.std(npcurv2)

				min2 = mean - std
				max2 = mean + std
		fig = plt.figure()
		xa = np.arange(-0.1, 0.1,0.001)
		plt.hist(npcurv1, bins=xa)
		plt.hist(npcurv2, bins=xa)
		plt.gca().set_yscale('log')
		plt.show()

		# Lookup table.
		lut = list()
		for idx in range(len(sources)):
			lut.append(vtk.vtkLookupTable())
			lut[idx].SetNumberOfColors(256)
			for i, item in enumerate(cc):
				lut[idx].SetTableValue(i, item[0], item[1], item[2], 1.0)
			if idx == 0:
				lut[idx].SetRange(min1, max1)
			if idx == 1:
				lut[idx].SetRange(min2, max2)
			if idx == 2:
				lut[idx].SetRange(min3, max3)
			if idx == 3:
				lut[idx].SetRange(np.min(parts), np.max(parts))
			lut[idx].Build()


		renderers = list()
		mappers = list()
		actors = list()
		textmappers = list()
		textactors = list()
		scalarbars = list()

		# Create a common text property.
		textProperty = vtk.vtkTextProperty()
		textProperty.SetFontSize(16)
		textProperty.SetJustificationToCentered()

		names = ['Sphere - Gaussian Curvature', 'Sphere - Mean Curvature', 'Sphere - Residuals to fitted sphere', 'Sphere - Residuals to fitted sphere']

		# Link the pipeline together. 
		for idx, item in enumerate(sources):
			#sources[idx].Update()


			mappers.append(vtk.vtkPolyDataMapper())
			mappers[idx].SetUseLookupTableScalarRange(1)
			if idx < 2:
				#print npcurv.min(), npcurv.max()
				#print npcurv

				mappers[idx].SetInputConnection(curvatures[idx].GetOutputPort())
				mappers[idx].SetLookupTable(lut[idx])
			else:
				#mappers[idx].SetInputConnection(sources[idx].GetOutputPort())
				mappers[idx].SetInputData(sources[idx])
				mappers[idx].SetColorModeToMapScalars()
				#mappers[idx].GetPointData().SetScalars(scalars)
				mappers[idx].SetLookupTable(lut[idx])

			actors.append(vtk.vtkActor())
			actors[idx].SetMapper(mappers[idx])

			textmappers.append(vtk.vtkTextMapper())
			textmappers[idx].SetInput(names[idx])
			textmappers[idx].SetTextProperty(textProperty)

			textactors.append(vtk.vtkActor2D())
			textactors[idx].SetMapper(textmappers[idx])
			textactors[idx].SetPosition(150, 16)

			scalarbars.append(vtk.vtkScalarBarActor())
			scalarbars[idx].SetLookupTable(lut[idx])

			renderers.append(vtk.vtkRenderer())

		gridDimensions = 2

		for idx in range(len(sources)):
			if idx < gridDimensions * gridDimensions:
				renderers.append(vtk.vtkRenderer)

		rendererSize = 300

		# Create the RenderWindow
		#
		renderWindow = vtk.vtkRenderWindow()
		renderWindow.SetSize(rendererSize * gridDimensions, 0)

		# Add and position the renders to the render window.
		viewport = list()
		for row in range(gridDimensions):
			for col in range(gridDimensions):
				idx = row * gridDimensions + col

				viewport[:] = []
				viewport.append(float(col) * rendererSize / (gridDimensions * rendererSize))
				viewport.append(float(gridDimensions - (row+1)) * rendererSize / (gridDimensions * rendererSize))
				viewport.append(float(col+1)*rendererSize / (gridDimensions * rendererSize))
				viewport.append(float(gridDimensions - row) * rendererSize / (gridDimensions * rendererSize))

				if idx > (len(sources) - 1):
					continue

				renderers[idx].SetViewport(viewport)
				renderWindow.AddRenderer(renderers[idx])

				renderers[idx].AddActor(actors[idx])
				renderers[idx].AddActor(textactors[idx])
				renderers[idx].AddActor(scalarbars[idx])
				renderers[idx].SetBackground(0.4,0.4,0.4)

		interactor = vtk.vtkRenderWindowInteractor()
		interactor.SetRenderWindow(renderWindow)

		renderWindow.Render()

		interactor.Start()

if __name__ == "__main__":
	po = CurvaturesDemo()
	po.CurvaturesDemo()

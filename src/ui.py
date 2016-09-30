import vtk
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

import numpy as np

class QVTKRenderWindowInteractorWheelfix(QVTKRenderWindowInteractor):
	def wheelEvent(self, ev):
		if ev.angleDelta().y() >= 0:
			self._Iren.MouseWheelForwardEvent()
		else:
			self._Iren.MouseWheelBackwardEvent()

class VTKMainWindow(QtWidgets.QMainWindow):
	def __init__(self, parent=None):
		QtWidgets.QMainWindow.__init__(self,parent)

		self.vtkCameras = []
	def addCamera(self, renderer, camera):
		vtkCamera = vtk.vtkCamera()
		camPosition = np.linalg.inv(camera.R).dot(camera.position)
		vtkCamera.SetPosition(camPosition[0], camPosition[1], camPosition[2])
		vtkCamera.SetFocalPoint(0,0,0)
		vtkCamera.SetThickness(500)
		planesArray = [0 for i in range(24)]

		vtkCamera.GetFrustumPlanes(camera.w/camera.h, planesArray)

		planes = vtk.vtkPlanes()
		planes.SetFrustumPlanes(planesArray)

		frustum = vtk.vtkFrustumSource()
		frustum.SetPlanes(planes)
		#cube.SetXLength(100)
		#cube.SetYLength(100)
		#cube.SetZLength(100)

		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputConnection(frustum.GetOutputPort())

		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		actor.GetProperty().SetRepresentationToWireframe()

		renderer.AddActor(actor)
		self.vtkCameras.append(vtkCamera)
	def setupFloorGrid(self, renderer, noCells, gridScale):
		grid = vtk.vtkRectilinearGrid()
 
		grid.SetDimensions(noCells[0],1,noCells[1]);
	 
		xArray = vtk.vtkDoubleArray()
		zArray = vtk.vtkDoubleArray()
		for i in range(noCells[0]):
			xArray.InsertNextValue(i*gridScale[0]);
		for i in range(noCells[1]):
			zArray.InsertNextValue(i*gridScale[1]);

		grid.SetXCoordinates(xArray);
		#grid.SetYCoordinates(yArray);
		grid.SetZCoordinates(zArray);

		#shrinkFilter = vtk.vtkShrinkFilter()
		#shrinkFilter.SetInputData(grid)

		rmapper = vtk.vtkDataSetMapper()
		#rmapper.SetInputConnection(shrinkFilter.GetOutputPort())
		rmapper.SetInputData(grid)

		ractor = vtk.vtkActor()
		ractor.GetProperty().SetRepresentationToWireframe()
		#ractor.SetOrigin(gridScale*gridSize[0], 0, gridScale*gridSize[2])
		ractor.SetPosition(-gridScale[0]*noCells[0]/2,-2000,-gridScale[1]*noCells[1]/2)
		ractor.GetProperty().SetColor(0.5,0.5,0.5)
		ractor.SetMapper(rmapper)

		renderer.AddActor(ractor)

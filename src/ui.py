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
	def __init__(self, stereoCameras, mainVTKRenderer, parent=None):
		QtWidgets.QMainWindow.__init__(self,parent)

		self.mainVTKRenderer = mainVTKRenderer

		self.cameraDock = QtWidgets.QDockWidget('Cameras', self)
		self.cameraDock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea);
		self.cameraList = QtWidgets.QListWidget(self.cameraDock)
		self.cameraDock.setWidget(self.cameraList)

		self.vtkCameras = []
		self.itemCameraMM = dict()
		for i,stereoCamera in stereoCameras.iteritems():
			self.addCamera(self.mainVTKRenderer,stereoCamera.A)
			self.addCamera(self.mainVTKRenderer,stereoCamera.B)
			self.addCamera(self.mainVTKRenderer,stereoCamera.C)

			listItem = QtWidgets.QListWidgetItem("Camera " + str(stereoCamera.name), self.cameraList)
			listItem.setFlags(listItem.flags() | QtCore.Qt.ItemIsUserCheckable)
			listItem.setCheckState(QtCore.Qt.Checked)
			self.cameraList.addItem(listItem)
			
			self.itemCameraMM[listItem] = stereoCamera

		self.cameraList.itemChanged.connect(self.onCameraToggle)

		self.stereoCameras = stereoCameras
		self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.cameraDock);

	def addCamera(self, renderer, camera):
		camPosition = np.linalg.inv(camera.R).dot(camera.position)
		vtkCamera = vtk.vtkCamera()
		vtkCamera.SetPosition(camPosition[0], camPosition[1], camPosition[2])
		vtkCamera.SetFocalPoint(0,0,0)
		vtkCamera.SetThickness(500)
		planesArray = [0 for i in range(24)]

		vtkCamera.GetFrustumPlanes(camera.w/camera.h, planesArray)

		planes = vtk.vtkPlanes()
		planes.SetFrustumPlanes(planesArray)

		frustum = vtk.vtkFrustumSource()
		frustum.SetPlanes(planes)
		'''
		cube = vtk.vtkCubeSource()
		cube.SetXLength(100)
		cube.SetYLength(100)
		cube.SetZLength(100)
		'''

		mapper = vtk.vtkPolyDataMapper()
		#mapper.SetInputConnection(cube.GetOutputPort())
		mapper.SetInputConnection(frustum.GetOutputPort())

		#transform = vtk.vtkPerspectiveTransform()
		#transform.SetupCamera(camPosition[0], camPosition[1], camPosition[2], 0, 0, 0, 0, 1, 0)

		actor = vtk.vtkActor()
		actor.SetMapper(mapper)

		'''
		print actor.GetXRange()
		print actor.GetYRange()
		print actor.GetZRange()
		'''
		#actor.SetUserMatrix(transform.GetMatrix())
		actor.SetPosition(camPosition[0], camPosition[1], camPosition[2])
		actor.GetProperty().SetRepresentationToWireframe()

		renderer.AddActor(actor)
		renderer.AddActor(labelActor)
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

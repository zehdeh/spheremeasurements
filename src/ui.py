import vtk
from math import sqrt
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import cv2

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
		self.settingsPanel = QtWidgets.QWidget()
		self.settingsDock = QtWidgets.QDockWidget('Properties',self)
		self.settingsDock.setWidget(self.settingsPanel)

		self.cameraActors = []
		self.setupSettings()

		self.vtkCameras = []
		self.itemCameraMM = dict()
		camPositions = vtk.vtkPoints()
		labels = vtk.vtkStringArray()
		labels.SetName('label')

		for i,stereoCamera in stereoCameras.iteritems():
			self.addCamera(self.mainVTKRenderer,stereoCamera.A)
			self.addCamera(self.mainVTKRenderer,stereoCamera.B)
			camPosition = self.addCamera(self.mainVTKRenderer,stereoCamera.C)

			labels.InsertNextValue(str(stereoCamera.name))
			camPositions.InsertNextPoint(camPosition[0], camPosition[1], camPosition[2])

			
		labelPolyData = vtk.vtkPolyData()
		labelPolyData.SetPoints(camPositions)
		labelPolyData.GetPointData().AddArray(labels)

		hier = vtk.vtkPointSetToLabelHierarchy()
		hier.SetInputData(labelPolyData)
		hier.SetLabelArrayName('label')

		labelMapper = vtk.vtkLabelPlacementMapper()
		labelMapper.SetInputConnection(hier.GetOutputPort())

		labelActor = vtk.vtkActor2D()
		labelActor.SetMapper(labelMapper)

		self.mainVTKRenderer.AddActor(labelActor)

		self.stereoCameras = stereoCameras
		self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.settingsDock)
	def setupSettings(self):
		grid = QtWidgets.QGridLayout()
		grid.setSpacing(10)

		displayMode = QtWidgets.QLabel('Display mode')
		errorMode = QtWidgets.QRadioButton('Voxels')
		errorMode.setChecked(True)
		errorMode.toggled.connect(lambda: self.switchDisplayMode(errorMode, 0))
		normalsMode = QtWidgets.QRadioButton('Normal vectors')
		normalsMode.toggled.connect(lambda: self.switchDisplayMode(normalsMode, 1))


		grid.addWidget(displayMode, 1, 0)
		grid.addWidget(errorMode, 1, 1)
		grid.addWidget(normalsMode, 1, 2)

		self.settingsPanel.setLayout(grid)

	def addCamera(self, renderer, camera):
		camPosition = np.linalg.inv(camera.R).dot(camera.position)

		vtkCamera = vtk.vtkCamera()
		vtkCamera.SetPosition(camPosition[0], camPosition[1], camPosition[2])

		rot, J = cv2.Rodrigues(camera.R)
		theta = sqrt(rot[0]**2 + rot[1]**2 + rot[2]**2)
		v = rot/theta

		vtkCamera.SetFocalPoint(0,0,0)
		#vtkCamera.SetDistance(camera.focalLength)
		vtkCamera.SetThickness(500)

		transform = vtk.vtkTransform()
		transform.RotateWXYZ(theta, v[0], v[1], v[2])
		#vtkCamera.SetUserTransform(transform)


		planesArray = [0 for i in range(24)]

		vtkCamera.GetFrustumPlanes(camera.w/camera.h, planesArray)
		#print 'planes:'
		#print np.min(planesArray)
		#print np.max(planesArray)

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
		actor.GetProperty().SetColor(0.4,0.4,0.4)
		actor.GetProperty().SetOpacity(0.5)
		#print 'Bounds:'
		#print np.min(actor.GetBounds())
		#print np.max(actor.GetBounds())

		'''
		print actor.GetXRange()
		print actor.GetYRange()
		print actor.GetZRange()
		'''
		#actor.SetUserMatrix(transform.GetMatrix())
		#actor.SetPosition(camPosition[0], camPosition[1], camPosition[2])
		actor.GetProperty().SetRepresentationToWireframe()

		renderer.AddActor(actor)
		self.cameraActors.append(actor)
		self.vtkCameras.append(vtkCamera)

		return camPosition
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

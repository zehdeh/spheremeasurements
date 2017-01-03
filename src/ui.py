import vtk
from math import sqrt
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
from vtk.util.numpy_support import numpy_to_vtk
from src.utils import rodrigues
import cv2

import numpy as np

class QVTKRenderWindowInteractorWheelfix(QVTKRenderWindowInteractor):
	def wheelEvent(self, ev):
		if ev.angleDelta().y() >= 0:
			self._Iren.MouseWheelForwardEvent()
		else:
			self._Iren.MouseWheelBackwardEvent()

class MainWindow(QtWidgets.QMainWindow):
	def __init__(self, stereoCameras, gridSize, gridScale, errorMatrix, vectorField, confidenceMatrix):
		QtWidgets.QMainWindow.__init__(self,parent=None)

		self._gridSize = gridSize
		self._gridScale = gridScale

		self._errorMatrix = errorMatrix
		self._vectorField = vectorField

		self._mainVTKRenderer = vtk.vtkRenderer()
		self._wireframeGridActor = None

		self.cameraActors = []
		self.setupSettings()

		self.vtkCameras = []
		self.itemCameraMM = dict()
		camPositions = vtk.vtkPoints()

		self._stereoCameras = stereoCameras

		tabWidget = QtWidgets.QTabWidget();
		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)

		self._interactor = self.errorGridViewer.GetRenderWindow().GetInteractor()
		interactorStyle = vtk.vtkInteractorStyleTerrain()
		self._interactor.SetInteractorStyle(interactorStyle)
		
		self.confidenceMatrix = confidenceMatrix
		self.errorGridViewer.GetRenderWindow().AddRenderer(self._mainVTKRenderer)
		self.setupVolume()
		self.setupErrorGrid()
		self.setupVectorField()
		tabWidget.addTab(self.errorGridViewer, 'Camera Visibility')

		self.setCentralWidget(tabWidget)
		self._interactor.Initialize()
		self.show()
	def setupCameras(self):
		labels = vtk.vtkStringArray()
		labels.SetName('label')

		for i,stereoCamera in self._stereoCameras.iteritems():
			self.addCamera(self._mainVTKRenderer,stereoCamera.A)
			self.addCamera(self._mainVTKRenderer,stereoCamera.B)
			camPosition = self.addCamera(self._mainVTKRenderer,stereoCamera.C)

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

		self._mainVTKRenderer.AddActor(labelActor)
	def setupVolume(self):

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0.0, 0.0)
		alphaChannelFunc.AddPoint(0.8, 1.0)

		colorFunc = vtk.vtkColorTransferFunction()
		colorFunc.AddRGBPoint(0.0, 1.0, 1.0, 1.0)
		colorFunc.AddRGBPoint(0.01, 0.706, 0.016, 0.150)

		scalarBar = vtk.vtkScalarBarActor()
		scalarBar.SetLookupTable(colorFunc)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.IndependentComponentsOff()
		volumeProperty.SetColor(0,colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)

		self.volumeMapper = vtk.vtkSmartVolumeMapper()

		self.volume = vtk.vtkVolume()
		self.volume.SetOrigin(self._gridSize[0]/2., self._gridSize[1]/2.,self._gridSize[2]/2.)
		self.volume.SetScale(self._gridScale[0], self._gridScale[1], self._gridScale[2])

		self.volume.SetPosition(0,0,0)
		self.volume.SetMapper(self.volumeMapper)
		self.volume.SetProperty(volumeProperty)


		self._mainVTKRenderer.AddVolume(self.volume)
		#self._mainVTKRenderer.AddActor(scalarBar)

	def onCameraToggle(self, item):
		self.itemCameraMM[item].visible = (item.checkState() == 2)
		self.rebuildCoverage()
		self.imageImport.Update()
		self.gaussianFilter.Update()
		self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())
		self.volumeMapper.Update()
		self.errorGridViewer.GetRenderWindow().Render()
		self._mainVTKRenderer.ResetCamera()
	def switchDisplayMode(self,btn,mode):
		if btn.isChecked():
			if mode == 0:
				self.volume.VisibilityOn()
				self.vectorActor.VisibilityOff()
			else:
				self.volume.VisibilityOff()
				self.vectorActor.VisibilityOn()
			self._mainVTKRenderer.GetRenderWindow().Render()
	def rebuildCoverage(self):
		alphaValues = np.zeros((self._gridSize[0], self._gridSize[1], self._gridSize[2]), dtype=np.float)
		colorValues = np.zeros((self._gridSize[0], self._gridSize[1], self._gridSize[2]), dtype=np.float)

		noCameras = float(len(self._stereoCameras))
		for i,stereoCamera in self._stereoCameras.iteritems():
			if stereoCamera.visible:
				pass
				#alphaValues += stereoCamera.visibilityMatrix / noCameras
		#alphaValues[:,:,:] = 1.0
		
		#colorValues[:,:,:] = -1.0
		
		colorValues = self._errorMatrix
			#colorValues = 255
			#alphaValues += 10*errorMatrix
		#colorValues = colorValues/np.max(colorValues)

		alphaValues = np.absolute(colorValues)
		#alphaValues = np.abs(colorValues)
		#alphaValues = self.confidenceMatrix

		alphaImporter = vtkImageImportFromArray()
		alphaImporter.SetArray(alphaValues)
		colorImporter = vtkImageImportFromArray()
		colorImporter.SetArray(colorValues)

		self.imageImport = vtk.vtkImageAppendComponents()
		self.imageImport.AddInputConnection(colorImporter.GetOutputPort())
		self.imageImport.AddInputConnection(alphaImporter.GetOutputPort())
	def setupVectorField(self):
		normals = self._vectorField.reshape(-1,3)
		nonzeroNormalIdx = np.unique(np.nonzero(normals)[0])
		normals = normals[nonzeroNormalIdx]

		indexPoints = np.ascontiguousarray(np.array(np.unravel_index(nonzeroNormalIdx, self._gridSize)).T)
		tmp = indexPoints[:,0].copy()
		indexPoints[:,0] = indexPoints[:,2]
		indexPoints[:,2] = tmp
		dataArray = numpy_to_vtk(indexPoints,array_type=vtk.VTK_INT)

		normalsArray = numpy_to_vtk(normals)

		points = vtk.vtkPoints()
		points.SetData(dataArray)
		polyData = vtk.vtkPolyData()
		polyData.SetPoints(points)
		polyData.GetPointData().SetNormals(normalsArray)

		arrowSource = vtk.vtkLineSource()

		glyph3D = vtk.vtkGlyph3D()
		glyph3D.SetSourceConnection(arrowSource.GetOutputPort())
		glyph3D.SetInputData(polyData)
		glyph3D.SetVectorModeToUseNormal()
		glyph3D.SetScaleModeToScaleByVector()
		glyph3D.SetScaleFactor(1)
		glyph3D.Update()

		vectorMapper = vtk.vtkPolyDataMapper()
		vectorMapper.SetInputConnection(glyph3D.GetOutputPort())

		self.vectorActor = vtk.vtkActor()
		self.vectorActor.VisibilityOff()
		self.vectorActor.SetMapper(vectorMapper)
		self.vectorActor.SetScale(self._gridScale[0], self._gridScale[1], self._gridScale[2])
		self.vectorActor.SetOrigin(self._gridSize[0]/2., self._gridSize[1]/2.,self._gridSize[2]/2.)

		self._mainVTKRenderer.AddActor(self.vectorActor)

	def setupErrorGrid(self):
			#errorGrid[0:(gridSize[0]-1), gridSize[1]-1-i,0:(gridSize[2]-1)] = i
		self.rebuildCoverage()
		self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())

	def setupSettings(self):
		self.settingsPanel = QtWidgets.QWidget()
		self.settingsDock = QtWidgets.QDockWidget('Properties',self)
		self.settingsDock.setWidget(self.settingsPanel)

		grid = QtWidgets.QGridLayout()
		grid.setSpacing(10)

		displayModeGroup = QtWidgets.QGroupBox('Display mode')
		errorMode = QtWidgets.QRadioButton('Voxels', displayModeGroup)
		errorMode.setChecked(True)
		errorMode.toggled.connect(lambda: self.switchDisplayMode(errorMode, 0))
		normalsMode = QtWidgets.QRadioButton('Normal vectors', displayModeGroup)
		normalsMode.toggled.connect(lambda: self.switchDisplayMode(normalsMode, 1))

		hbox = QtWidgets.QHBoxLayout()
		hbox.addWidget(errorMode)
		hbox.addWidget(normalsMode)
		displayModeGroup.setLayout(hbox)

		grid.addWidget(displayModeGroup, 1, 0)

		wireFrameGridGroup = QtWidgets.QGroupBox('Show wireframe grid')
		noWfButton = QtWidgets.QRadioButton('None', wireFrameGridGroup)
		noWfButton.setChecked(True)
		noWfButton.toggled.connect(lambda(active): active and self.toggleWireFrameGrid(0))
		floorWfButton = QtWidgets.QRadioButton('Floor only', wireFrameGridGroup)
		floorWfButton.toggled.connect(lambda(active): active and self.toggleWireFrameGrid(1))
		fullWfButton = QtWidgets.QRadioButton('Full', wireFrameGridGroup)
		fullWfButton.toggled.connect(lambda(active): active and self.toggleWireFrameGrid(2))

		hbox = QtWidgets.QHBoxLayout()
		hbox.addWidget(noWfButton)
		hbox.addWidget(floorWfButton)
		hbox.addWidget(fullWfButton)
		wireFrameGridGroup.setLayout(hbox)

		grid.addWidget(wireFrameGridGroup, 2, 0)

		self.settingsPanel.setLayout(grid)
		self.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.settingsDock)

	def addCamera(self, renderer, camera):
		camPosition = np.linalg.inv(camera.R).dot(camera.position)

		vtkCamera = vtk.vtkCamera()
		vtkCamera.SetPosition(camPosition[0], camPosition[1], camPosition[2])

		#print 'Rotation matrix:'
		#print camera.R
		#print 'Comparison:'
		rot, J = cv2.Rodrigues(camera.R)
		#print rot
		#print rodrigues(camera.R)
		#rot = rodrigues(camera.R)
		theta = sqrt(rot[0]**2 + rot[1]**2 + rot[2]**2)
		v = rot/theta

		vtkCamera.SetFocalPoint(0,0,0)
		#vtkCamera.SetDistance(camera.focalLength)

		transform = vtk.vtkTransform()
		transform.RotateWXYZ(theta, v[0], v[1], v[2])
		#vtkCamera.SetUserTransform(transform)

		vtkCamera.SetThickness(500)


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

		self._mainVTKRenderer.AddActor(actor)
		self.cameraActors.append(actor)
		self.vtkCameras.append(vtkCamera)

		return camPosition
	def toggleWireFrameGrid(self, mode):
		if self._wireframeGridActor is not None:
			self._mainVTKRenderer.RemoveActor(self._wireframeGridActor)
		
		if mode == 1:
			self.setupWireframeGrid(True)
		elif mode == 2:
			self.setupWireframeGrid(False)
		self.errorGridViewer.GetRenderWindow().Render()

	def setupWireframeGrid(self, floorGrid=False):
		grid = vtk.vtkRectilinearGrid()
 
		if floorGrid:
			grid.SetDimensions(self._gridSize[0],1,self._gridSize[2]);
	 
			xArray = numpy_to_vtk(np.arange(self._gridSize[0])*self._gridScale[0], array_type=vtk.VTK_INT)
			zArray = numpy_to_vtk(np.arange(self._gridSize[2])*self._gridScale[2], array_type=vtk.VTK_INT)

			grid.SetXCoordinates(xArray);
			grid.SetZCoordinates(zArray);
		else:
			grid.SetDimensions(self._gridSize[0],self._gridSize[1],self._gridSize[2]);

			xArray = numpy_to_vtk(np.arange(self._gridSize[0])*self._gridScale[0], array_type=vtk.VTK_INT)
			yArray = numpy_to_vtk(np.arange(self._gridSize[1])*self._gridScale[1], array_type=vtk.VTK_INT)
			zArray = numpy_to_vtk(np.arange(self._gridSize[2])*self._gridScale[2], array_type=vtk.VTK_INT)

			grid.SetXCoordinates(xArray);
			grid.SetYCoordinates(yArray);
			grid.SetZCoordinates(zArray);

		rmapper = vtk.vtkDataSetMapper()
		rmapper.SetInputData(grid)

		self._wireframeGridActor = vtk.vtkActor()
		self._wireframeGridActor.GetProperty().SetRepresentationToWireframe()
		if floorGrid:
			self._wireframeGridActor.SetPosition(-self._gridScale[0]*self._gridSize[0]/2,-(self._gridScale[1]*self._gridSize[1])/2,-self._gridScale[2]*self._gridSize[1]/2)
		else:
			self._wireframeGridActor.SetPosition(-self._gridScale[0]*self._gridSize[0]/2,-self._gridScale[1]*self._gridSize[1]/2,-self._gridScale[2]*self._gridSize[1]/2)
		self._wireframeGridActor.GetProperty().SetOpacity(0.5)
		self._wireframeGridActor.GetProperty().SetColor(0.5,0.5,0.5)
		self._wireframeGridActor.SetMapper(rmapper)

		self._mainVTKRenderer.AddActor(self._wireframeGridActor)

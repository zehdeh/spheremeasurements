import vtk
from math import sqrt
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from scipy.ndimage.filters import convolve
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
		self._confidenceMatrix = confidenceMatrix

		self._mainVTKRenderer = vtk.vtkRenderer()
		self._wireframeGridActor = None

		self.cameraActors = []
		self.setupSettings()

		self.vtkCameras = []
		self.itemCameraMM = dict()

		self._stereoCameras = stereoCameras

		tabWidget = QtWidgets.QTabWidget();
		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)

		self._interactor = self.errorGridViewer.GetRenderWindow().GetInteractor()
		interactorStyle = vtk.vtkInteractorStyleTerrain()
		self._interactor.SetInteractorStyle(interactorStyle)
		
		self.errorGridViewer.GetRenderWindow().AddRenderer(self._mainVTKRenderer)
		self.setupVolume()
		self.setupErrorGrid()
		self.setupVectorField()
		tabWidget.addTab(self.errorGridViewer, 'Camera Visibility')

		axes = vtk.vtkAxesActor()
		axes.SetTotalLength(1000,1000,1000)
		self._mainVTKRenderer.AddActor(axes)

		self.setCentralWidget(tabWidget)
		self._interactor.Initialize()
		self.show()

		self.setupCameras()
	def setupCameras(self):
		labels = vtk.vtkStringArray()
		labels.SetName('label')

		camPositions = vtk.vtkPoints()
		for i,stereoCamera in self._stereoCameras.iteritems():
			print stereoCamera.name
			self.addCamera(self._mainVTKRenderer,stereoCamera.A)
			self.addCamera(self._mainVTKRenderer,stereoCamera.B)
			self.addCamera(self._mainVTKRenderer,stereoCamera.C)
			camPosition = np.linalg.inv(stereoCamera.C.R).dot(stereoCamera.C.position)

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

		#self._mainVTKRenderer.AddActor(labelActor)
	def setupVolume(self):

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0.0, 0.0)
		alphaChannelFunc.AddPoint(1.0, 0.7)

		mean = self._errorMatrix[np.nonzero(self._errorMatrix)].mean()
		std = self._errorMatrix[np.nonzero(self._errorMatrix)].std()
		upperBound = mean + std-0.02
		lowerBound = 0
		print upperBound

		colorFunc = vtk.vtkColorTransferFunction()
		colorFunc.AddRGBPoint(lowerBound, 0.0, 0.0, 1.0)
		colorFunc.AddRGBPoint((upperBound - lowerBound)/2, 0.0, 1.0, 0.0)
		colorFunc.AddRGBPoint(upperBound, 1.0, 0.0, 0.0)

		scalarBar = vtk.vtkScalarBarActor()
		scalarBar.SetLookupTable(colorFunc)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.IndependentComponentsOff()
		volumeProperty.SetColor(colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)
		volumeProperty.ShadeOff()

		self.volumeMapper = vtk.vtkSmartVolumeMapper()

		self.volume = vtk.vtkVolume()
		self.volume.SetOrigin(self._gridSize[0]/2., self._gridSize[1]/2.,self._gridSize[2]/2.)
		self.volume.SetScale(self._gridScale[0], self._gridScale[1], self._gridScale[2])

		self.volume.SetPosition(0,0,0)
		self.volume.SetMapper(self.volumeMapper)
		self.volume.SetProperty(volumeProperty)

		self._mainVTKRenderer.AddVolume(self.volume)
		#self._mainVTKRenderer.AddActor(scalarBar)

	def switchDisplayMode(self,btn,mode):
		if btn.isChecked():
			if mode == 0:
				self.volume.VisibilityOn()
				self.vectorActor.VisibilityOff()
			else:
				self.volume.VisibilityOff()
				self.vectorActor.VisibilityOn()
			self._mainVTKRenderer.GetRenderWindow().Render()
	def rebuildCoverage(self, convolutionPasses = 0, growConvolution = True):
		colorValues = self._errorMatrix
		alphaValues = self._confidenceMatrix

		for i in range(convolutionPasses):
			sizeX = 3+i if growConvolution else 3
			boxFilter = np.zeros((sizeX,sizeX,sizeX))
			boxFilter[:,:,:] = (1./float(sizeX**3))

			colorValues = convolve(colorValues, boxFilter)
			alphaValues = convolve(alphaValues, boxFilter)

			colorValues[np.nonzero(self._errorMatrix)] = self._errorMatrix[np.nonzero(self._errorMatrix)]

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
		self.rebuildCoverage(2,True)
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

		convolutionGroup = QtWidgets.QGroupBox('Convolution')
		convolutionPassesLabel = QtWidgets.QLabel('Passes:', convolutionGroup)
		convolutionPassesSpinner = QtWidgets.QSpinBox(convolutionGroup)
		convolutionPassesSpinner.setMinimum(0)
		convolutionPassesSpinner.setValue(2)
		convolutionGrowButton = QtWidgets.QCheckBox('Grow', convolutionGroup)
		convolutionGrowButton.setCheckState(2)
		convolutionApplyButton = QtWidgets.QPushButton('Apply', convolutionGroup)

		def applyConvolution(passes, grow):
			self.rebuildCoverage(passes,grow)
			self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())
			self.volumeMapper.Update()
			self.errorGridViewer.GetRenderWindow().Render()
		convolutionApplyButton.clicked.connect(lambda: applyConvolution(convolutionPassesSpinner.value(), convolutionGrowButton.checkState() == 2))

		vbox = QtWidgets.QHBoxLayout()
		vbox.addWidget(convolutionPassesLabel)
		vbox.addWidget(convolutionPassesSpinner)
		vbox.addWidget(convolutionGrowButton)
		vbox.addWidget(convolutionApplyButton)
		convolutionGroup.setLayout(vbox)

		grid.addWidget(convolutionGroup, 3, 0)

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

		planes = vtk.vtkPlanes()
		planes.SetFrustumPlanes(planesArray)

		frustum = vtk.vtkFrustumSource()
		frustum.SetPlanes(planes)
		#print frustum.GetOutput().GetBounds()
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
		if -3000 > np.min(actor.GetBounds()) or 3000 < np.max(actor.GetBounds()):
			print 'Warning: camera dropped because it was out of bounds'
			print frustum.GetOutput().GetBounds()
			return

		#actor.SetUserMatrix(transform.GetMatrix())
		#actor.SetPosition(camPosition[0], camPosition[1], camPosition[2])
		actor.GetProperty().SetRepresentationToWireframe()

		self._mainVTKRenderer.AddActor(actor)
		self.cameraActors.append(actor)
		self.vtkCameras.append(vtkCamera)

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

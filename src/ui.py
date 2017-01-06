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

VIEWMODE = type('enum', (), dict(error=0,curvature=1))
class MainWindow(QtWidgets.QMainWindow):
	def __init__(self, stereoCameras, gridSize, gridScale, errorMatrix, curvatureMatrix, vectorField, confidenceMatrix, verbose=False):
		QtWidgets.QMainWindow.__init__(self,parent=None)

		self._gridSize = gridSize
		self._gridScale = gridScale

		self._errorMatrix = errorMatrix
		self._curvatureMatrix = curvatureMatrix
		self._vectorField = vectorField
		self._confidenceMatrix = confidenceMatrix

		self._verbose = verbose

		self._mainVTKRenderer = vtk.vtkRenderer()
		self._wireframeGridActor = None

		self.cameraActors = []
		self.setupSettings()

		self.vtkCameras = []
		self.itemCameraMM = dict()

		self._stereoCameras = stereoCameras

		self._viewMode = VIEWMODE.error

		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(self)

		self._interactor = self.errorGridViewer.GetRenderWindow().GetInteractor()
		interactorStyle = vtk.vtkInteractorStyleTerrain()
		self._interactor.SetInteractorStyle(interactorStyle)
		
		self.errorGridViewer.GetRenderWindow().AddRenderer(self._mainVTKRenderer)
		self.setupVolume()
		self.setupErrorGrid()
		self.setupVectorField()

		self.viewDock = QtWidgets.QDockWidget('View',self)
		self.viewChooserToolBox = QtWidgets.QToolBox()
		self.analysisListWidget = QtWidgets.QListWidget(self.viewChooserToolBox)

		errorViewButton = QtWidgets.QListWidgetItem('Error')
		curvatureViewButton = QtWidgets.QListWidgetItem('Curvature')

		def callback(e):
			if e == errorViewButton:
				self._viewMode = VIEWMODE.error

				self.repaintVolume()
			else:
				self._viewMode = VIEWMODE.curvature

				self.repaintVolume()

		self.analysisListWidget.addItem(errorViewButton)
		self.analysisListWidget.addItem(curvatureViewButton)
		self.analysisListWidget.itemClicked.connect(lambda e: callback(e))

		self.analysisListWidget.setCurrentItem(errorViewButton)

		self.cameraListWidget = QtWidgets.QListWidget(self.viewChooserToolBox)

		self.modelListWidget = QtWidgets.QListWidget(self.viewChooserToolBox)
		self.modelListWidget.addItem('Curvature Model')

		self.viewChooserToolBox.addItem(self.analysisListWidget, 'Analysis')
		self.viewChooserToolBox.addItem(self.cameraListWidget, 'Camera Coverage')
		self.viewChooserToolBox.addItem(self.modelListWidget, 'Model')


		self.viewDock.setWidget(self.viewChooserToolBox)
		self.viewDock.setFeatures(QtWidgets.QDockWidget.NoDockWidgetFeatures)

		self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.viewDock)

		axes = vtk.vtkAxesActor()
		axes.SetTotalLength(1000,1000,1000)

		axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(20)
		axes.GetXAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
		axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(20)
		axes.GetYAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
		axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(20)
		axes.GetZAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()

		self._mainVTKRenderer.AddActor(axes)

		self.setCentralWidget(self.errorGridViewer)
		self._interactor.Initialize()

		dw = QtWidgets.QDesktopWidget()

		self.resize(dw.width()*0.7, dw.height()*0.7)
		self.move((dw.width()-dw.width()*0.7)/2, (dw.height()-dw.height()*0.7)/2)

		self.show()

		self.setupCameras()

	def setupColorFunction(self):
		plottedMatrix = self._errorMatrix if self._viewMode == VIEWMODE.error else self._curvatureMatrix

		mean = self._errorMatrix[np.nonzero(plottedMatrix)].mean()
		std = self._errorMatrix[np.nonzero(plottedMatrix)].std()/2.5
		upperBound = mean + std
		lowerBound = 0
		if self._verbose:
			print 'Color range (' + str(lowerBound) +',' + str(upperBound) + ')'

		self._colorFunc = vtk.vtkColorTransferFunction()
		self._colorFunc.AddRGBPoint(lowerBound, 0.0, 0.0, 1.0)
		self._colorFunc.AddRGBPoint((upperBound - lowerBound)/2, 1.0, 1.0, 1.0)
		self._colorFunc.AddRGBPoint(upperBound, 1.0, 0.0, 0.0)

		self._scalarBar.SetLookupTable(self._colorFunc)

	def setupCameras(self):
		labels = vtk.vtkStringArray()
		labels.SetName('label')

		camPositions = vtk.vtkPoints()
		for i,stereoCamera in self._stereoCameras.iteritems():
			self.addCamera(self._mainVTKRenderer,stereoCamera.A)
			self.addCamera(self._mainVTKRenderer,stereoCamera.B)
			self.addCamera(self._mainVTKRenderer,stereoCamera.C)
			camPosition = np.linalg.inv(stereoCamera.C.R).dot(stereoCamera.C.position)

			labels.InsertNextValue(str(stereoCamera.name))
			camPositions.InsertNextPoint(camPosition[0], camPosition[1], camPosition[2])

			camListItem = QtWidgets.QListWidgetItem('Camera ' + str(stereoCamera.name))
			camListItem.setFlags(camListItem.flags() | QtCore.Qt.ItemIsUserCheckable);
			camListItem.setCheckState(QtCore.Qt.Checked)

			self.cameraListWidget.addItem(camListItem)

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
		alphaChannelFunc.AddPoint(1.0, 0.7)

		self._scalarBar = vtk.vtkScalarBarActor()
		self._scalarBar.GetLabelTextProperty().ItalicOff()
		self._scalarBar.GetLabelTextProperty().BoldOff()
		self._scalarBar.SetNumberOfLabels(3)
		self._scalarBar.AnnotationTextScalingOff()
		self._scalarBar.SetMaximumWidthInPixels(100)
		self._scalarBar.SetPosition(0.9, 0.15)
		self._mainVTKRenderer.AddActor(self._scalarBar)

		self.setupColorFunction()

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.IndependentComponentsOff()
		volumeProperty.SetColor(self._colorFunc)
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

	def switchDisplayMode(self,btn,mode):
		if btn.isChecked():
			if mode == 0:
				self.volume.VisibilityOn()
				self.vectorActor.VisibilityOff()
			else:
				self.volume.VisibilityOff()
				self.vectorActor.VisibilityOn()
			self._mainVTKRenderer.GetRenderWindow().Render()
	def repaintVolume(self):
		self.generateVolumePlot()

		self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())
		self.volumeMapper.Update()

		self.setupColorFunction()

		self.errorGridViewer.GetRenderWindow().Render()
	def generateVolumePlot(self, convolutionPasses = 0, growConvolution = True):
		colorValues = self._errorMatrix if self._viewMode == VIEWMODE.error else self._curvatureMatrix
		alphaValues = self._confidenceMatrix

		growConvolution = self._convolutionGrowButton.checkState() == QtCore.Qt.Checked
		convolutionPasses = self._convolutionPassesSpinner.value()

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

		arrowSource = vtk.vtkArrowSource()
		arrowSource.SetShaftResolution(3)
		arrowSource.SetTipResolution(3)

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
		self.vectorActor.GetProperty().SetRepresentationToWireframe()

		self._mainVTKRenderer.AddActor(self.vectorActor)

	def setupErrorGrid(self):
		self.generateVolumePlot()
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

		convolutionGroup = QtWidgets.QGroupBox('Box Filter')
		convolutionPassesLabel = QtWidgets.QLabel('Passes:', convolutionGroup)
		self._convolutionPassesSpinner = QtWidgets.QSpinBox(convolutionGroup)
		self._convolutionPassesSpinner.setMinimum(0)
		self._convolutionPassesSpinner.setValue(2)
		self._convolutionGrowButton = QtWidgets.QCheckBox('Grow', convolutionGroup)
		self._convolutionGrowButton.setCheckState(2)
		convolutionApplyButton = QtWidgets.QPushButton('Apply', convolutionGroup)

		convolutionApplyButton.clicked.connect(lambda: self.repaintVolume())

		vbox = QtWidgets.QHBoxLayout()
		vbox.addWidget(convolutionPassesLabel)
		vbox.addWidget(self._convolutionPassesSpinner)
		vbox.addWidget(self._convolutionGrowButton)
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

		outOfBounds = True
		while outOfBounds:
			frustum = vtk.vtkFrustumSource()
			frustum.SetPlanes(planes)

			mapper = vtk.vtkPolyDataMapper()
			mapper.SetInputConnection(frustum.GetOutputPort())

			#transform = vtk.vtkPerspectiveTransform()
			#transform.SetupCamera(camPosition[0], camPosition[1], camPosition[2], 0, 0, 0, 0, 1, 0)

			actor = vtk.vtkActor()
			actor.SetMapper(mapper)
			actor.GetProperty().SetColor(0.4,0.4,0.4)
			actor.GetProperty().SetOpacity(0.5)
			outOfBounds = -3000 > np.min(actor.GetBounds()) or 3000 < np.max(actor.GetBounds())

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

		xArray = numpy_to_vtk(np.arange(self._gridSize[0])*self._gridScale[0], array_type=vtk.VTK_INT)
		zArray = numpy_to_vtk(np.arange(self._gridSize[2])*self._gridScale[2], array_type=vtk.VTK_INT)
		if floorGrid:
			grid.SetDimensions(self._gridSize[0],1,self._gridSize[2]);
	 

			grid.SetXCoordinates(xArray);
			grid.SetZCoordinates(zArray);

		else:
			yArray = numpy_to_vtk(np.arange(self._gridSize[1])*self._gridScale[1], array_type=vtk.VTK_INT)

			grid.SetDimensions(self._gridSize[0],self._gridSize[1],self._gridSize[2]);

			grid.SetXCoordinates(xArray);
			grid.SetYCoordinates(yArray);
			grid.SetZCoordinates(zArray);

		mapper = vtk.vtkDataSetMapper()
		mapper.SetInputData(grid)

		self._wireframeGridActor = vtk.vtkActor()
		self._wireframeGridActor.GetProperty().SetRepresentationToWireframe()
		if floorGrid:
			self._wireframeGridActor.SetPosition(-self._gridScale[0]*self._gridSize[0]/2,-(self._gridScale[1]*self._gridSize[1])/2,-self._gridScale[2]*self._gridSize[1]/2)
		else:
			self._wireframeGridActor.SetPosition(-self._gridScale[0]*self._gridSize[0]/2,-self._gridScale[1]*self._gridSize[1]/2,-self._gridScale[2]*self._gridSize[1]/2)
		self._wireframeGridActor.GetProperty().SetOpacity(0.5)
		self._wireframeGridActor.GetProperty().SetColor(0.5,0.5,0.5)
		self._wireframeGridActor.SetMapper(mapper)

		self._mainVTKRenderer.AddActor(self._wireframeGridActor)

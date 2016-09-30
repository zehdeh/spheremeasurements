#! /usr/bin/python

import sys
from math import atan2,sqrt
import os

import numpy as np

import cv2
import vtk
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from src.OBJIO import loadOBJ, writeOBJ
from src.fitting import fitSphere
from src.thirdparty.body.loaders.scanner_frame import Pod
from src.mesh import Mesh
from src.ui import VTKMainWindow, QVTKRenderWindowInteractorWheelfix
from opendr.camera import ProjectPoints

class StereoCamera(object):
	def __init__(self):
		self.A = None
		self.B = None
		self.C = None
		self.ppointsA = None
		self.ppointsB = None
		self.name = None
		self.visibilityMatrix = None
		self.visible = True

class MainWindow(VTKMainWindow):
	def __init__(self, stereoCameras, centerPoints, gridSize, gridScale, parent = None):
		VTKMainWindow.__init__(self, parent)

		tabWidget = QtWidgets.QTabWidget();
		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)
		
		cameraDock = QtWidgets.QDockWidget('Cameras', self)
		cameraDock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea);
		
		cameraList = QtWidgets.QListWidget()
		cameraDock.setWidget(cameraList)

		self.stereoCameras = stereoCameras
		self.setupErrorGrid(self.errorGridViewer,stereoCameras, centerPoints, gridSize, gridScale, cameraList)
		tabWidget.addTab(self.errorGridViewer, 'Camera Visibility')

		cameraList.itemChanged.connect(self.onCameraToggle)

		self.addDockWidget(QtCore.Qt.RightDockWidgetArea, cameraDock);
		self.setCentralWidget(tabWidget)
		self.show()
	
	def onCameraToggle(self, item):
		self.itemCameraMM[item].visible = (item.checkState() == 2)
		self.rebuildCoverage()
		self.volumeMapper.Update()
		self.errorGridViewer.GetRenderWindow().Render()
	def rebuildCoverage(self):
		errorGrid = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)

		for i,stereoCamera in self.stereoCameras.iteritems():
			if stereoCamera.visible:
				errorGrid += stereoCamera.visibilityMatrix
			else:
				print 'not using this one'
		#errorGrid *= 100

		dataString = errorGrid.tostring()
		self.gridDataImporter.CopyImportVoidPointer(dataString, len(dataString))
	def setupErrorGrid(self,errorGridViewer,stereoCameras, centerPoints, gridSize, gridScale, uiCameraList):
		self.gridRenderer = vtk.vtkRenderer()
		errorGridViewer.GetRenderWindow().AddRenderer(self.gridRenderer)
		iren = errorGridViewer.GetRenderWindow().GetInteractor()

		self.setupFloorGrid(self.gridRenderer, [gridSize[0], gridSize[2]], [gridScale[0], gridScale[2]])


		self.itemCameraMM = dict()
		for i,stereoCamera in stereoCameras.iteritems():
			self.addCamera(self.gridRenderer,stereoCamera.A)
			self.addCamera(self.gridRenderer,stereoCamera.B)
			self.addCamera(self.gridRenderer,stereoCamera.C)

			listItem = QtWidgets.QListWidgetItem("Camera " + str(stereoCamera.name), uiCameraList)
			listItem.setFlags(listItem.flags() | QtCore.Qt.ItemIsUserCheckable)
			listItem.setCheckState(QtCore.Qt.Checked)
			uiCameraList.addItem(listItem)
			
			self.itemCameraMM[listItem] = stereoCamera
		#print errorGrid
			
			#errorGrid[0:(gridSize[0]-1), gridSize[1]-1-i,0:(gridSize[2]-1)] = i
		#for cp in centerPoints:
		#	cp[0] += gridSize[0]*gridScale/2
		#	cp[1] += gridSize[1]*gridScale/2
		#	cp[2] += gridSize[2]*gridScale/2
		#	i = int(cp[0] / gridScale)
		#	j = int(cp[1] / gridScale)
		#	k = int(cp[2] / gridScale)
		#	errorGrid[i,j,k] += 150

		self.gridDataImporter = vtk.vtkImageImport()
		self.rebuildCoverage()

		self.gridDataImporter.SetDataScalarTypeToUnsignedChar()
		self.gridDataImporter.SetNumberOfScalarComponents(1)

		self.gridDataImporter.SetDataExtent(0, (gridSize[0] - 1), 0, (gridSize[1] - 1), 0, (gridSize[2] - 1))
		self.gridDataImporter.SetWholeExtent(0, (gridSize[0] - 1), 0, (gridSize[1] - 1), 0, (gridSize[2] - 1))
		self.gridDataImporter.SetDataSpacing(1.0,1.0,1.0)
		self.gridDataImporter.SetDataOrigin(0,0,0)

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0, 0.0)
		alphaChannelFunc.AddPoint(255, 0.2)

		colorFunc = vtk.vtkColorTransferFunction()
		colorFunc.AddRGBPoint(0, 0,0,1.0)
		colorFunc.AddRGBPoint(50, 0,0,1.0)
		colorFunc.AddRGBPoint(255, 0.0,0,1.0)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.SetColor(colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)

		self.volumeMapper = vtk.vtkSmartVolumeMapper()
		self.volumeMapper.SetInputConnection(self.gridDataImporter.GetOutputPort())

		volume = vtk.vtkVolume()
		volume.SetOrigin(gridSize[0]/2, gridSize[1]/2,gridSize[2]/2)
		volume.SetScale(gridScale[0], gridScale[1], gridScale[2])
		#volume.SetPosition(0,-gridSize[1]*gridScale/2,0)
		#volume.SetPosition(0,0,0)
		volume.SetMapper(self.volumeMapper)
		volume.SetProperty(volumeProperty)

		#interactorStyle = vtk.vtkInteractorStyleTerrain()
		#iren.SetInteractorStyle(interactorStyle)
		self.gridRenderer.AddVolume(volume)
		iren.Initialize()

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('2 Arguments needed')
		sys.exit(0)
	app = QtWidgets.QApplication(sys.argv)

	folderPath = sys.argv[1]
	files = os.listdir(folderPath)
	centerPoints = []
	for fileName in files:
		if fileName.endswith('.obj'):
			vertices, faces, normals = loadOBJ(folderPath + '/' + fileName)
			bounds = Mesh(vertices, faces, normals).getBounds()
			p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
			cp, radius = fitSphere(vertices,p0,150, bounds)
			centerPoints.append(cp)

	calibrationFolderPath = sys.argv[2]
	files = os.listdir(calibrationFolderPath)
	cameras = []
	i = 0
	for fileName in files:
		if fileName.endswith('.tka'):
			pod = Pod(calibrationFolderPath + fileName,'',image_scale=1.0,to_meters=False,bg_image_filename='')
			cameras.append((fileName[0:-4],pod.camera()))
			i += 1

	gridSize = [50,50,50]
	scannerVolumeSize = [3000,3000,3000]
	gridScale = [scannerVolumeSize[0] / gridSize[0], scannerVolumeSize[1] / gridSize[1], scannerVolumeSize[2] / gridSize[2]]

	cameras = [(cam[0],cam[1]) for cam in cameras]

	stereoCameras = dict()
	#stereoCameras[list(headIndices)] = []
	
	for cam in cameras:
		camNo = int(cam[0][0:2])
		camSign = cam[0][3]
		if camNo not in stereoCameras:
			stereoCameras[camNo] = StereoCamera()
			stereoCameras[camNo].name = camNo

		if camSign == 'A':
			stereoCameras[camNo].A = cam[1]
			stereoCameras[camNo].ppointsA = ProjectPoints(f=cam[1].f.ravel(), rt=cam[1].r.ravel(), t=cam[1].t.ravel(), k=cam[1].k.ravel(), c=cam[1].c.ravel())
		elif camSign == 'B':
			stereoCameras[camNo].B = cam[1]
			stereoCameras[camNo].ppointsB = ProjectPoints(f=cam[1].f.ravel(), rt=cam[1].r.ravel(), t=cam[1].t.ravel(), k=cam[1].k.ravel(), c=cam[1].c.ravel())
		else:
			stereoCameras[camNo].C = cam[1]

	for l,stereoCamera in stereoCameras.iteritems():
		fileName = calibrationFolderPath + '/coverageMatrices/' + str(stereoCamera.name) + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2])
		if os.path.isfile(fileName + '.npy'):
			stereoCamera.visibilityMatrix = np.load(fileName + '.npy')
		else:
			visibilityMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)
			print 'processing camera...'
			for i in range(gridSize[0]):
				for j in range(gridSize[1]):
					for k in range(gridSize[2]):
						x = (i*gridScale[0]) - (gridSize[0]*gridScale[0]/2)
						y = (j*gridScale[1]) - (gridSize[1]*gridScale[1]/2)
						z = (k*gridScale[2]) - (gridSize[2]*gridScale[2]/2)

						stereoCamera.ppointsA.v = [x,y,z]
						stereoCamera.ppointsB.v = [x,y,z]

						x1 = stereoCamera.ppointsA.r[0]
						y1 = stereoCamera.ppointsA.r[1]
						x2 = stereoCamera.ppointsB.r[0]
						y2 = stereoCamera.ppointsB.r[1]

						if x1 < 1600 and x1 >= 0 and x2 < 1600 and x2 > 0 and y1 < 1200 and y1 >= 0 and y2 < 1200 and y2 >= 0:
							visibilityMatrix[k,j,i] = 1
			stereoCamera.visibilityMatrix = visibilityMatrix
			if not os.path.exists(calibrationFolderPath + '/coverageMatrices/'):
				os.makedirs(calibrationFolderPath + '/coverageMatrices/')
			np.save(fileName, visibilityMatrix)
	
	window = MainWindow(stereoCameras, centerPoints, gridSize, gridScale)
	sys.exit(app.exec_())

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
from src.fitting import fitSphere, fittingErrorSphere
from src.thirdparty.body.loaders.scanner_frame import Pod
from src.mesh import Mesh
from src.ui import VTKMainWindow, QVTKRenderWindowInteractorWheelfix
from opendr.camera import ProjectPoints
from src.calibration import getStereoCamerasFromCalibration, StereoCamera


class MainWindow(VTKMainWindow):
	def __init__(self, stereoCameras, centerPoints, rmses, gridSize, gridScale, parent = None):
		VTKMainWindow.__init__(self, stereoCameras, vtk.vtkRenderer(), parent)

		tabWidget = QtWidgets.QTabWidget();
		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)
		
		
		self.errorGridViewer.GetRenderWindow().AddRenderer(self.mainVTKRenderer)
		self.setupErrorGrid(stereoCameras, self.mainVTKRenderer, centerPoints, gridSize, gridScale)
		tabWidget.addTab(self.errorGridViewer, 'Camera Visibility')


		self.setCentralWidget(tabWidget)
		self.show()
	
	def onCameraToggle(self, item):
		self.itemCameraMM[item].visible = (item.checkState() == 2)
		self.rebuildCoverage()
		self.volumeMapper.Update()
		self.errorGridViewer.GetRenderWindow().Render()
	def rebuildCoverage(self):
		errorGrid = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=(np.uint8,2))
		#colorValues = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)
		for i in range(gridSize[1]):
			errorGrid[:,i,:,0] = 0
			errorGrid[:,i,:,1] = 0

		for i,stereoCamera in self.stereoCameras.iteritems():
			if stereoCamera.visible:
				pass
				#alphaValues += stereoCamera.visibilityMatrix
		#alphaValues *= 10

		dataString = errorGrid.tostring()
		self.gridDataImporter.CopyImportVoidPointer(dataString, len(dataString))
	def setupErrorGrid(self,stereoCameras, renderer, centerPoints, gridSize, gridScale):
		iren = self.errorGridViewer.GetRenderWindow().GetInteractor()

		self.setupFloorGrid(self.mainVTKRenderer, [gridSize[0], gridSize[2]], [gridScale[0], gridScale[2]])

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
		self.gridDataImporter.SetNumberOfScalarComponents(2)

		self.gridDataImporter.SetDataExtent(0, (gridSize[0] - 1), 0, (gridSize[1] - 1), 0, (gridSize[2] - 1))
		self.gridDataImporter.SetWholeExtent(0, (gridSize[0] - 1), 0, (gridSize[1] - 1), 0, (gridSize[2] - 1))
		self.gridDataImporter.SetDataSpacing(1.0,1.0,1.0)
		self.gridDataImporter.SetDataOrigin(0,0,0)

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0, 1.0)
		alphaChannelFunc.AddPoint(255, 0.0)

		colorFunc = vtk.vtkColorTransferFunction()
		colorFunc.AddRGBPoint(0, 1.0,1.0,1.0)
		colorFunc.AddRGBPoint(50, 0.5,0.5,0.5)
		colorFunc.AddRGBPoint(255, 1.0,1.0,1.0)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.SetColor(0,colorFunc)
		volumeProperty.SetScalarOpacity(1,alphaChannelFunc)

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
		self.mainVTKRenderer.AddVolume(volume)
		iren.Initialize()

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('2 Arguments needed')
		sys.exit(0)
	app = QtWidgets.QApplication(sys.argv)

	folderPath = sys.argv[1]
	files = os.listdir(folderPath)
	centerPoints = []
	rmses = []
	for fileName in files:
		if fileName.endswith('.obj'):
			vertices, faces, normals = loadOBJ(folderPath + '/' + fileName)
			bounds = Mesh(vertices.T, faces, normals).getBounds()
			p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
			cp, radius = fitSphere(vertices,cp,150, bounds)
			centerPoints.append(cp)
			rmses.append(np.linalg.norm(fittingErrorSphere(p0,vertices))/np.sqrt(len(vertices)))

	gridSize = [50,50,50]
	scannerVolumeSize = [3000,3000,3000]
	gridScale = [scannerVolumeSize[0] / gridSize[0], scannerVolumeSize[1] / gridSize[1], scannerVolumeSize[2] / gridSize[2]]
	
	calibrationFolderPath = sys.argv[2]
	stereoCameras = getStereoCamerasFromCalibration(calibrationFolderPath)

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
	
	window = MainWindow(stereoCameras, centerPoints, rmses, gridSize, gridScale)
	sys.exit(app.exec_())

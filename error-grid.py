#! /usr/bin/python

import timeit
import sys
from math import atan2,sqrt,fabs
import os

import numpy as np

import cv2
import vtk
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtk.util.vtkImageImportFromArray import vtkImageImportFromArray
from src.OBJIO import loadOBJ, writeOBJ, loadOBJviaVTK
from src.fitting import fitSphere, fittingErrorSphere, calculateMeanCurvature
from src.thirdparty.body.loaders.scanner_frame import Pod
from src.mesh import Mesh
from src.ui import VTKMainWindow, QVTKRenderWindowInteractorWheelfix
from opendr.camera import ProjectPoints
from src.calibration import getStereoCamerasFromCalibration, StereoCamera


class MainWindow(VTKMainWindow):
	def __init__(self, stereoCameras, centerPoints, errorMatrix, gridSize, gridScale, parent = None):
		VTKMainWindow.__init__(self, stereoCameras, vtk.vtkRenderer(), parent)

		tabWidget = QtWidgets.QTabWidget();
		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)
		
		self.centerPoints = centerPoints
		self.errorMatrix = errorMatrix
		self.errorGridViewer.GetRenderWindow().AddRenderer(self.mainVTKRenderer)
		self.setupErrorGrid(stereoCameras, self.mainVTKRenderer, gridSize, gridScale)
		tabWidget.addTab(self.errorGridViewer, 'Camera Visibility')


		self.setCentralWidget(tabWidget)
		self.show()
	
	def onCameraToggle(self, item):
		self.itemCameraMM[item].visible = (item.checkState() == 2)
		self.rebuildCoverage()
		self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())
		self.volumeMapper.Update()
		self.errorGridViewer.GetRenderWindow().Render()
		self.mainVTKRenderer.ResetCamera()
	def rebuildCoverage(self):
		alphaValues = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
		colorValues = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)

		noCameras = float(len(self.stereoCameras))
		for i,stereoCamera in self.stereoCameras.iteritems():
			if stereoCamera.visible:
				pass
				#alphaValues += stereoCamera.visibilityMatrix / noCameras
		#alphaValues[:,:,:] = 1.0
		
		#colorValues[:,:,:] = -1.0
		
		colorValues = self.errorMatrix
			#colorValues = 255
			#alphaValues += 10*errorMatrix
		colorValues = colorValues/np.max(colorValues)

		alphaValues = np.abs(colorValues)
		#alphaValues = np.fabs(colorValues)
		'''
		for cp in self.centerPoints:
			x = cp[0] + gridSize[0]*gridScale[0]/2
			y = cp[1] + gridSize[1]*gridScale[1]/2
			z = cp[2] + gridSize[2]*gridScale[2]/2
			i = int(x / gridScale[0])
			j = int(y / gridScale[1])
			k = int(z / gridScale[2])
			colorValues[i,j,k] = 255#min(255,colorValues[i,j,k] + 150)
			alphaValues[i,j,k] += 50
		'''

		alphaImporter = vtkImageImportFromArray()
		alphaImporter.SetArray(alphaValues)
		colorImporter = vtkImageImportFromArray()
		colorImporter.SetArray(colorValues)

		self.imageImport = vtk.vtkImageAppendComponents()
		self.imageImport.AddInputConnection(colorImporter.GetOutputPort())
		self.imageImport.AddInputConnection(alphaImporter.GetOutputPort())

	def setupErrorGrid(self,stereoCameras, renderer, gridSize, gridScale):
		iren = self.errorGridViewer.GetRenderWindow().GetInteractor()

		self.setupFloorGrid(self.mainVTKRenderer, [gridSize[0], gridSize[2]], [gridScale[0], gridScale[2]])

		#print errorGrid
			
			#errorGrid[0:(gridSize[0]-1), gridSize[1]-1-i,0:(gridSize[2]-1)] = i

		self.rebuildCoverage()

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0, 0.0)
		alphaChannelFunc.AddPoint(1.0, 0.8)

		colorFunc = vtk.vtkColorTransferFunction()
		#colorFunc.AddRGBPoint(-1.0, 0.230, 0.299, 0.754)
		colorFunc.AddRGBPoint(0.0, 1.0, 1.0, 1.0)
		colorFunc.AddRGBPoint(1.0, 0.706, 0.016, 0.150)

		scalarBar = vtk.vtkScalarBarActor()
		scalarBar.SetLookupTable(colorFunc)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.IndependentComponentsOff()
		volumeProperty.SetColor(0,colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)

		self.volumeMapper = vtk.vtkSmartVolumeMapper()
		self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())

		volume = vtk.vtkVolume()
		volume.SetOrigin(gridSize[0]/2, gridSize[1]/2,gridSize[2]/2)
		volume.SetScale(gridScale[0], gridScale[1], gridScale[2])
		#volume.SetPosition(0,-gridSize[1]*gridScale/2,0)
		#volume.SetPosition(0,0,0)
		volume.SetMapper(self.volumeMapper)
		volume.SetProperty(volumeProperty)

		interactorStyle = vtk.vtkInteractorStyleTerrain()
		iren.SetInteractorStyle(interactorStyle)
		self.mainVTKRenderer.AddVolume(volume)
		self.mainVTKRenderer.AddActor(scalarBar)
		iren.Initialize()

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('2 Arguments needed')
		sys.exit(0)
	app = QtWidgets.QApplication(sys.argv)

	folderPath = sys.argv[1]
	files = os.listdir(folderPath)
	centerPoints = []

	#gridSize = [50,50,50]
	gridSize = [100,100,100]
	scannerVolumeSize = [3000,3000,3000]
	gridScale = [scannerVolumeSize[0] / gridSize[0], scannerVolumeSize[1] / gridSize[1], scannerVolumeSize[2] / gridSize[2]]
	print files[230]

	matrixFileName = folderPath + '/meanCurvature' + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2]) + '.npy'
	if os.path.isfile(matrixFileName + '.npy'):
		totalErrorMatrix = np.load(matrixFileName + '.npy')
	else:
		errorMatrices = []
		for fileName in files:
			if fileName.endswith('.obj'):
				print 'Processing ' + fileName[0:-4]
				vertices, faces, normals, polyData = loadOBJviaVTK(folderPath + fileName)

				#print len(faces), len(faces2)
				bounds = Mesh(vertices.T, faces, normals).getBounds()
				p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
				cp, radius = fitSphere(vertices,p0,150, bounds)
				centerPoints.append(cp)

				#errors = fittingErrorSphere(cp.tolist()+[radius],vertices) - radius
				errors = calculateMeanCurvature(polyData)

				errorMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
				for v,e in zip(vertices, errors):
					x = v[0] + gridSize[0]*gridScale[0]/2
					y = v[1] + gridSize[1]*gridScale[1]/2
					z = v[2] + gridSize[2]*gridScale[2]/2
					i = int(x / gridScale[0])
					j = int(y / gridScale[1])
					k = int(z / gridScale[2])
					errorMatrix[i,j,k] = e
				errorMatrix = np.fabs(errorMatrix)
				#print errorMatrix.max(), errorMatrix.min()
				errorMatrices.append(errorMatrix)
		
		totalErrorMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
		for errorMatrix in errorMatrices:
			totalErrorMatrix += errorMatrix
		np.save(matrixFileName, totalErrorMatrix)

	calibrationFolderPath = sys.argv[2]
	stereoCameras = getStereoCamerasFromCalibration(calibrationFolderPath)

	for l,stereoCamera in stereoCameras.iteritems():
		stereoCamera.visibilityMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)
		'''
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
			'''
	
	window = MainWindow(stereoCameras, centerPoints, totalErrorMatrix, gridSize, gridScale)
	sys.exit(app.exec_())

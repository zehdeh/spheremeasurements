#! /usr/bin/python

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
import pymesh


class MainWindow(VTKMainWindow):
	def __init__(self, stereoCameras, centerPoints, errorMatrix, vectorField, gridSize, gridScale, parent = None):
		VTKMainWindow.__init__(self, stereoCameras, vtk.vtkRenderer(), parent)

		tabWidget = QtWidgets.QTabWidget();
		self.errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)
		
		self.centerPoints = centerPoints
		self.errorMatrix = errorMatrix
		self.vectorField = vectorField
		self.errorGridViewer.GetRenderWindow().AddRenderer(self.mainVTKRenderer)
		self.setupErrorGrid(stereoCameras, self.mainVTKRenderer, gridSize, gridScale)
		tabWidget.addTab(self.errorGridViewer, 'Camera Visibility')


		self.setCentralWidget(tabWidget)
		self.show()
	
	def onCameraToggle(self, item):
		self.itemCameraMM[item].visible = (item.checkState() == 2)
		self.rebuildCoverage()
		self.imageImport.Update()
		self.gaussianFilter.Update()
		self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())
		self.volumeMapper.Update()
		self.errorGridViewer.GetRenderWindow().Render()
		self.mainVTKRenderer.ResetCamera()
	def switchDisplayMode(self,btn,mode):
		if btn.isChecked():
			if mode == 0:
				self.volume.VisibilityOn()
			else:
				self.volume.VisibilityOff()
			self.mainVTKRenderer.GetRenderWindow().Render()
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

		#self.setupFloorGrid(self.mainVTKRenderer, [gridSize[0], gridSize[2]], [gridScale[0], gridScale[2]])

		#print errorGrid
			
			#errorGrid[0:(gridSize[0]-1), gridSize[1]-1-i,0:(gridSize[2]-1)] = i
		normals = self.vectorField.reshape(-1,3)

		normalsArray = vtk.vtkDoubleArray()
		normalsArray.SetNumberOfComponents(3)

		arrowSource = vtk.vtkArrowSource()
		glyph3D = vtk.vtkGlyph3D()

		self.rebuildCoverage()

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0.003, 0.0)
		alphaChannelFunc.AddPoint(0.11, 0.8)

		colorFunc = vtk.vtkColorTransferFunction()
		#colorFunc.AddRGBPoint(-1.0, 0.230, 0.299, 0.754)
		colorFunc.AddRGBPoint(0.0, 1.0, 1.0, 1.0)
		colorFunc.AddRGBPoint(0.01, 0.706, 0.016, 0.150)

		scalarBar = vtk.vtkScalarBarActor()
		scalarBar.SetLookupTable(colorFunc)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.IndependentComponentsOff()
		volumeProperty.SetColor(0,colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)

		self.gaussianFilter = vtk.vtkImageGaussianSmooth()
		self.gaussianFilter.SetInputConnection(self.imageImport.GetOutputPort())
		self.gaussianFilter.Update()

		self.volumeMapper = vtk.vtkSmartVolumeMapper()
		#self.volumeMapper.SetInputConnection(self.imageImport.GetOutputPort())
		self.volumeMapper.SetInputConnection(self.gaussianFilter.GetOutputPort())

		self.volume = vtk.vtkVolume()
		self.volume.SetOrigin(gridSize[0]/2., gridSize[1]/2.,gridSize[2]/2.)
		self.volume.SetScale(gridScale[0], gridScale[1], gridScale[2])
		#volume.SetPosition(0,-gridSize[1]*gridScale/2,0)
		self.volume.SetPosition(0,0,0)
		self.volume.SetMapper(self.volumeMapper)
		self.volume.SetProperty(volumeProperty)

		cubeAxesActor = vtk.vtkCubeAxesActor()
		cubeAxesActor.SetBounds(self.volume.GetBounds())

		cubeAxesActor.SetFlyMode(vtk.VTK_FLY_FURTHEST_TRIAD)
		cubeAxesActor.SetGridLineLocation(vtk.VTK_GRID_LINES_FURTHEST)
		cubeAxesActor.DrawXGridlinesOn()
		cubeAxesActor.DrawYGridlinesOn()
		cubeAxesActor.DrawZGridlinesOn()

		cubeAxesActor.XAxisMinorTickVisibilityOff()
		cubeAxesActor.YAxisMinorTickVisibilityOff()
		cubeAxesActor.ZAxisMinorTickVisibilityOff()

		cubeAxesActor.XAxisLabelVisibilityOff()
		cubeAxesActor.YAxisLabelVisibilityOff()
		cubeAxesActor.ZAxisLabelVisibilityOff()

		cubeAxesActor.SetCamera(self.mainVTKRenderer.GetActiveCamera())

		interactorStyle = vtk.vtkInteractorStyleTerrain()
		iren.SetInteractorStyle(interactorStyle)
		self.mainVTKRenderer.AddVolume(self.volume)
		#self.mainVTKRenderer.AddActor(cubeAxesActor)
		self.mainVTKRenderer.AddActor(scalarBar)
		iren.Initialize()
		
		self.mainVTKRenderer.ResetCamera()
		#self.mainVTKRenderer.GetActiveCamera().SetFocalPoint(0,0,0)
		#self.mainVTKRenderer.GetActiveCamera().SetPosition(0,0,-3000)

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('2 Arguments needed')
		sys.exit(0)
	app = QtWidgets.QApplication(sys.argv)

	folderPath = sys.argv[1]
	files = os.listdir(folderPath)
	centerPoints = []

	#gridSize = [50,50,50]
	#gridSize = [100,100,100]
	gridSize = [200,200,200]
	scannerVolumeSize = [3000,3000,3000]
	gridScale = [scannerVolumeSize[0] / gridSize[0], scannerVolumeSize[1] / gridSize[1], scannerVolumeSize[2] / gridSize[2]]

	matrixFileName = folderPath + '/meanCurvature' + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2])
	#matrixFileName = folderPath + '/fittingError' + '_' + str(gridSize[0]) + '_' + str(gridSize[1]) + '_' + str(gridSize[2])
	if os.path.isfile(matrixFileName + '.npy'):
		totalErrorMatrix = np.load(matrixFileName + '.npy')
		totalVectorField = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=(np.float,3))
	else:
		totalErrorMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
		totalVectorField = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=(np.float,3))
		nMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.int32)
		for fileName in files:
			if fileName.endswith('.obj'):
				print 'Processing ' + fileName[0:-4]
				vertices, faces, normals, polyData = loadOBJviaVTK(os.path.join(folderPath,fileName))

				#print len(faces), len(faces2)
				bounds = Mesh(vertices.T, faces, normals).getBounds()
				p0 = [bounds[0][0],bounds[1][0],bounds[2][0],150]
				cp, radius = fitSphere(vertices,p0,150, bounds)
				centerPoints.append(cp)

				#errors = fittingErrorSphere(cp.tolist()+[radius],vertices) - radius
				#mesh = pymesh.form_mesh(vertices,faces)
				#mesh.add_attribute('vertex_mean_curvature')
				#errors = mesh.get_attribute('vertex_mean_curvature')
				errors = calculateMeanCurvature(polyData)

				errorMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.float)
				vectorField = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=(np.float,3))
				for v,e,n in zip(vertices, errors, normals):
					x = v[0] + gridSize[0]*gridScale[0]/2
					y = v[1] + gridSize[1]*gridScale[1]/2
					z = v[2] + gridSize[2]*gridScale[2]/2
					k = int(x / gridScale[0])
					j = int(y / gridScale[1])
					i = int(z / gridScale[2])
					errorMatrix[i,j,k] = e
					vectorField[i,j,k] += e*n
					nMatrix[i,j,k] += 1
				errorMatrix = np.fabs(errorMatrix)
				#print errorMatrix.max(), errorMatrix.min()
				totalErrorMatrix += errorMatrix
				totalVectorField += vectorField
		
		totalErrorMatrix[np.nonzero(nMatrix)] = totalErrorMatrix[np.nonzero(nMatrix)] / nMatrix[np.nonzero(nMatrix)]
		np.save(matrixFileName, totalErrorMatrix)
		#print totalVectorField[:,:,:].shape

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
	
	window = MainWindow(stereoCameras, centerPoints, totalErrorMatrix, totalVectorField, gridSize, gridScale)
	sys.exit(app.exec_())

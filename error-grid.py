#! /usr/bin/python3

import sys
import os
import cv2

import numpy as np

import vtk
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from src.OBJIO import loadOBJ, writeOBJ
from src.fitting import fitSphere
from src.thirdparty.body.loaders.scanner_frame import Pod
from src.mesh import Mesh

class QVTKRenderWindowInteractorWheelfix(QVTKRenderWindowInteractor):
	def wheelEvent(self, ev):
		if ev.angleDelta().y() >= 0:
			self._Iren.MouseWheelForwardEvent()
		else:
			self._Iren.MouseWheelBackwardEvent()

class MainWindow(QtWidgets.QMainWindow):
	def __init__(self, cameras, centerPoints, parent = None):
		QtWidgets.QMainWindow.__init__(self,parent)

		tabWidget = QtWidgets.QTabWidget();
		errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)
		
		self.setupGrid(errorGridViewer,cameras, centerPoints)
		tabWidget.addTab(errorGridViewer, 'Camera Visibility')

		self.setCentralWidget(tabWidget)
		self.show()
	def setupGrid(self,errorGridViewer,cameras, centerPoints):
		renderer = vtk.vtkRenderer()
		errorGridViewer.GetRenderWindow().AddRenderer(renderer)
		iren = errorGridViewer.GetRenderWindow().GetInteractor()

		for camera in cameras:
			cube = vtk.vtkCubeSource()
			cube.SetXLength(100)
			cube.SetYLength(100)
			cube.SetZLength(100)

			mapper = vtk.vtkPolyDataMapper()
			mapper.SetInputConnection(cube.GetOutputPort())

			actor = vtk.vtkActor()
			actor.SetMapper(mapper)
			
			camPosition = -np.linalg.inv(camera.R).dot(camera.position)

			actor.SetPosition(camPosition[0], camPosition[1], camPosition[2])
			realR, J = cv2.Rodrigues(camera.R)
			realR = realR*180/np.pi
			actor.SetOrientation(realR[0], realR[1], realR[2])

			renderer.AddActor(actor)

		gridSize = 20
		gridScale = 30

		errorGrid = np.zeros([gridSize, gridSize, gridSize], dtype=np.uint8)
		for i in range(gridSize):
			errorGrid[0:(gridSize-1), gridSize-1-i,0:(gridSize-1)] = i
		#for cp in centerPoints:
		#	cp += cp[0] + gridSize*gridScale/2
		#	i = cp[0] / gridScale
		#	j = cp[1] / gridScale
		#	k = cp[2] / gridScale
		#	errorGrid[i,j,k] = 200

		dataImporter = vtk.vtkImageImport()
		dataString = errorGrid.tostring()
		dataImporter.CopyImportVoidPointer(dataString, len(dataString))
		dataImporter.SetDataScalarTypeToUnsignedChar()
		dataImporter.SetNumberOfScalarComponents(1)

		dataImporter.SetDataExtent(0, (gridSize - 1), 0, (gridSize - 1), 0, (gridSize - 1))
		dataImporter.SetWholeExtent(0, (gridSize - 1), 0, (gridSize - 1), 0, (gridSize - 1))

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0, 0)
		alphaChannelFunc.AddPoint(gridSize - 1, 0.1)

		colorFunc = vtk.vtkColorTransferFunction()
		colorFunc.AddRGBPoint(0, 0,0,1.0)
		colorFunc.AddRGBPoint(50, 0,0,1.0)
		colorFunc.AddRGBPoint(255, 0.0,0,1.0)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.SetColor(colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)

		volumeMapper = vtk.vtkSmartVolumeMapper()
		volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

		volume = vtk.vtkVolume()
		volume.SetOrigin(gridSize/2, gridSize/2,gridSize/2)
		volume.SetScale(gridScale)
		volume.SetPosition(0,-200,0)
		print(volume.GetBounds())
		volume.SetMapper(volumeMapper)
		volume.SetProperty(volumeProperty)

		renderer.AddVolume(volume)
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

	folderPath = sys.argv[2]
	files = os.listdir(folderPath)
	cameras = []
	i = 0
	for fileName in files:
		if fileName.endswith('.tka'):
			pod = Pod(folderPath + fileName,'',image_scale=1.0,to_meters=False,bg_image_filename='')
			cameras.append(pod.camera())
			i += 1

	window = MainWindow(cameras, centerPoints)
	sys.exit(app.exec_())

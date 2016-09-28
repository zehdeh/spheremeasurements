#! /usr/bin/python

import sys
from math import atan2,sqrt
import os

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
	
	def addCamera(self, renderer, camera):
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
		camera.lookat(np.array([0,0,0]))
		#rx = atan2(camera.R[2,1],camera.R[2,2])
		#ry = atan2(-camera.R[2,0],sqrt(camera.R[2,1]**2 + camera.R[2,2]**2))
		#rz = atan2(camera.R[1,0], camera.R[0,0])
		#print(rx)
		#realR = realR*180/np.pi
		#actor.SetOrientation(rx, ry, rz)

		renderer.AddActor(actor)
	def setupGrid(self,errorGridViewer,cameras, centerPoints):
		renderer = vtk.vtkRenderer()
		errorGridViewer.GetRenderWindow().AddRenderer(renderer)
		iren = errorGridViewer.GetRenderWindow().GetInteractor()

		for camera in cameras:
			self.addCamera(renderer,camera)

		gridSize = [50,50,50]
		gridScale = 60

		errorGrid = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)
		#for i in range(gridSize[1]):
		#	errorGrid[0:(gridSize[0]-1), gridSize[1]-1-i,0:(gridSize[2]-1)] = i
		for cp in centerPoints:
			cp[0] += gridSize[0]*gridScale/2
			cp[1] += gridSize[1]*gridScale/2
			cp[2] += gridSize[2]*gridScale/2
			i = int(cp[0] / gridScale)
			j = int(cp[1] / gridScale)
			k = int(cp[2] / gridScale)
			errorGrid[i,j,k] += 150

		dataImporter = vtk.vtkImageImport()
		dataString = errorGrid.tostring()
		dataImporter.CopyImportVoidPointer(dataString, len(dataString))
		dataImporter.SetDataScalarTypeToUnsignedChar()
		dataImporter.SetNumberOfScalarComponents(1)

		dataImporter.SetWholeExtent(0, (gridSize[0] - 1), 0, (gridSize[1] - 1), 0, (gridSize[2] - 1))
		dataImporter.SetDataExtentToWholeExtent()

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0, 0)
		alphaChannelFunc.AddPoint(gridSize[1] - 1, 0.5)

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
		volume.SetOrigin(gridSize[0]/2, gridSize[1]/2,gridSize[2]/2)
		volume.SetScale(gridScale)
		#volume.SetPosition(0,-gridSize[1]*gridScale/2,0)
		volume.SetPosition(0,0,0)
		volume.SetMapper(volumeMapper)
		volume.SetProperty(volumeProperty)

		#interactorStyle = vtk.vtkInteractorStyleTerrain()
		#iren.SetInteractorStyle(interactorStyle)
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

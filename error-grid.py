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
from src.ui import VTKMainWindow
from opendr.camera import ProjectPoints

class QVTKRenderWindowInteractorWheelfix(QVTKRenderWindowInteractor):
	def wheelEvent(self, ev):
		if ev.angleDelta().y() >= 0:
			self._Iren.MouseWheelForwardEvent()
		else:
			self._Iren.MouseWheelBackwardEvent()

class MainWindow(VTKMainWindow):
	def __init__(self, cameras, centerPoints, parent = None):
		VTKMainWindow.__init__(self, parent)

		tabWidget = QtWidgets.QTabWidget();
		errorGridViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)
		
		self.setupErrorGrid(errorGridViewer,cameras, centerPoints)
		tabWidget.addTab(errorGridViewer, 'Camera Visibility')

		self.setCentralWidget(tabWidget)
		self.show()
	
	def setupErrorGrid(self,errorGridViewer,cameras, centerPoints):
		renderer = vtk.vtkRenderer()
		errorGridViewer.GetRenderWindow().AddRenderer(renderer)
		iren = errorGridViewer.GetRenderWindow().GetInteractor()

		for camera in cameras:
			self.addCamera(renderer,camera[1])

		gridSize = [50,50,50]
		gridScale = 60

		self.setupFloorGrid(renderer, gridSize[0], gridScale)
		bwCameras = [(int(cam[0][0:2]),cam[1]) for cam in cameras if cam[0][-1] in ['A', 'B']]

		headIndices = set([cam[0] for cam in bwCameras])

		stereoCameras = dict.fromkeys(list(headIndices))
		#stereoCameras[list(headIndices)] = []
		
		for cam in bwCameras:
			ppoints = ProjectPoints(f=cam[1].f.ravel(), rt=cam[1].r.ravel(), t=cam[1].t.ravel(), k=cam[1].k.ravel(), c=cam[1].c.ravel())
			if stereoCameras[cam[0]] is None:
				stereoCameras[cam[0]] = [(cam[1],ppoints)]
			else:
				stereoCameras[cam[0]].append((cam[1],ppoints))

		errorGrid = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)
		for i in range(gridSize[0]):
			print i
			for j in range(gridSize[1]):
				for k in range(gridSize[2]):
					x = (i*gridScale) - (gridSize[0]*gridScale/2)
					y = (j*gridScale) - (gridSize[0]*gridScale/2)
					z = (k*gridScale) - (gridSize[0]*gridScale/2)

					for l,stereoCamera in stereoCameras.iteritems():
						stereoCamera[0][1].v = [x,y,z]
						stereoCamera[1][1].v = [x,y,z]

						x1 = stereoCamera[0][1].r[0]
						y1 = stereoCamera[0][1].r[1]
						x2 = stereoCamera[1][1].r[0]
						y2 = stereoCamera[1][1].r[1]

						if x1 < 1600 and x1 >= 0 and x2 < 1600 and x2 > 0 and y1 < 1200 and y1 >= 0 and y2 < 1200 and y2 >= 0:
							errorGrid[i,j,k] += 10
			#errorGrid[0:(gridSize[0]-1), gridSize[1]-1-i,0:(gridSize[2]-1)] = i
		#for cp in centerPoints:
		#	cp[0] += gridSize[0]*gridScale/2
		#	cp[1] += gridSize[1]*gridScale/2
		#	cp[2] += gridSize[2]*gridScale/2
		#	i = int(cp[0] / gridScale)
		#	j = int(cp[1] / gridScale)
		#	k = int(cp[2] / gridScale)
		#	errorGrid[i,j,k] += 150

		dataImporter = vtk.vtkImageImport()
		dataString = errorGrid.tostring()
		dataImporter.CopyImportVoidPointer(dataString, len(dataString))
		dataImporter.SetDataScalarTypeToUnsignedChar()
		dataImporter.SetNumberOfScalarComponents(1)

		dataImporter.SetWholeExtent(0, (gridSize[0] - 1), 0, (gridSize[1] - 1), 0, (gridSize[2] - 1))
		dataImporter.SetDataExtentToWholeExtent()

		alphaChannelFunc = vtk.vtkPiecewiseFunction()
		alphaChannelFunc.AddPoint(0, 0)
		alphaChannelFunc.AddPoint(255, 0.5)

		colorFunc = vtk.vtkColorTransferFunction()
		colorFunc.AddRGBPoint(0, 0,0,1.0)
		colorFunc.AddRGBPoint(50, 0,0,1.0)
		colorFunc.AddRGBPoint(255, 0.0,0,1.0)

		volumeProperty = vtk.vtkVolumeProperty()
		volumeProperty.SetColor(colorFunc)
		volumeProperty.SetScalarOpacity(alphaChannelFunc)

		volumeMapper = vtk.vtkSmartVolumeMapper()
		volumeMapper.SetInputConnection(dataImporter.GetOutputPort())
		#volumeMapper.SetInputDataObject(structuredGrid)

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
			cameras.append((fileName[0:-4],pod.camera()))
			i += 1

	window = MainWindow(cameras, centerPoints)
	sys.exit(app.exec_())

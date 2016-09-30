#! /usr/bin/python

import os
import sys
import vtk
import numpy as np
import cv2
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from src.thirdparty.body.loaders.scanner_frame import Pod
from opendr.camera import ProjectPoints
from src.ui import VTKMainWindow, QVTKRenderWindowInteractorWheelfix

from src.OBJIO import loadOBJ, writeOBJ

class MainWindow(VTKMainWindow):
	def __init__(self, cameras, vertices, faces, parent = None):
		VTKMainWindow.__init__(self,parent)

		tabWidget = QtWidgets.QTabWidget();
		errorTracerViewer = QVTKRenderWindowInteractorWheelfix(tabWidget)

		self.setupErrorTracer(errorTracerViewer,cameras, vertices, faces)
		tabWidget.addTab(errorTracerViewer, 'Camera Visibility')

		self.setCentralWidget(tabWidget)
		self.show()
	def setupErrorTracer(self, errorTracerViewer, cameras, vertices, faces):
		renderer = vtk.vtkRenderer()
		errorTracerViewer.GetRenderWindow().AddRenderer(renderer)
		iren = errorTracerViewer.GetRenderWindow().GetInteractor()

		for camera in cameras:
			print camera[0]
			self.addCamera(renderer,camera[1])

		self.setupFloorGrid(renderer, [50,50], [60,60])
		points = vtk.vtkPoints()
		for v in vertices:
			points.InsertNextPoint(v[0], v[1], v[2])

		polygons = vtk.vtkCellArray()
		for f in faces:
			polygon = vtk.vtkPolygon()
			polygon.GetPointIds().SetNumberOfIds(3)
			polygon.GetPointIds().SetId(0, f[0])
			polygon.GetPointIds().SetId(1, f[1])
			polygon.GetPointIds().SetId(2, f[2])

			polygons.InsertNextCell(polygon)

		colors = vtk.vtkUnsignedCharArray()
		colors.SetNumberOfComponents(3)
		for v in vertices:
			camera = cameras[7][1]
			ppoints = ProjectPoints(f=camera.f.ravel(), rt=camera.r.ravel(), t=camera.t.ravel(), k=camera.k.ravel(), c=camera.c.ravel())
			ppoints.v = v
			if ppoints.r[0] < 1600 and ppoints.r[0] >= 0 and ppoints.r[1] < 1200 and ppoints.r[1] >= 0:
				c = [255,0,0]
			else:
				c = [255,255,255]
			colors.InsertNextTupleValue(c)

		polydata = vtk.vtkPolyData()
		polydata.SetPoints(points)
		polydata.SetPolys(polygons)
		polydata.GetPointData().SetScalars(colors)

		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputData(polydata)

		actor = vtk.vtkActor()
		actor.SetMapper(mapper)

		renderer.AddActor(actor)

		iren.Initialize()

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print('2 Arguments needed')
		sys.exit(0)
	app = QtWidgets.QApplication(sys.argv)
	vertices, faces, normals = loadOBJ(sys.argv[1])

	folderPath = sys.argv[2]
	files = os.listdir(folderPath)
	cameras = []
	for fileName in files:
		if fileName.endswith('.tka'):
			imagePath = folderPath + fileName[0:-4] + '.bmp'
			pod = Pod(folderPath + fileName, imagePath,image_scale=1.0,to_meters=False,bg_image_filename=imagePath)
			cameras.append((fileName[0:-4],pod.camera()))

	window = MainWindow(cameras, vertices, faces)
	sys.exit(app.exec_())

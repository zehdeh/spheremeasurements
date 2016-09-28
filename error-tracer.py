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

from src.OBJIO import loadOBJ, writeOBJ

class QVTKRenderWindowInteractorWheelfix(QVTKRenderWindowInteractor):
	def wheelEvent(self, ev):
		if ev.angleDelta().y() >= 0:
			self._Iren.MouseWheelForwardEvent()
		else:
			self._Iren.MouseWheelBackwardEvent()

class MainWindow(QtWidgets.QMainWindow):
	def __init__(self, cameras, vertices, faces, parent = None):
		QtWidgets.QMainWindow.__init__(self,parent)

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
			self.addCamera(renderer,camera)

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
			ppoints = ProjectPoints(f=cameras[0].f.ravel(), rt=cameras[0].r.ravel(), t=cameras[0].t.ravel(), k=cameras[0].k.ravel(), c=cameras[0].c.ravel())
			ppoints.v = v
			print ppoints.r
			c = [255,0,0]
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

	def addCamera(self, renderer, camera):
		cube = vtk.vtkCubeSource()
		cube.SetXLength(100)
		cube.SetYLength(100)
		cube.SetZLength(100)

		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputConnection(cube.GetOutputPort())

		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		
		camPosition = np.linalg.inv(camera.R).dot(camera.position)

		actor.SetPosition(camPosition[0], camPosition[1], camPosition[2])
		camera.lookat(np.array([0,0,0]))
		#rx = atan2(camera.R[2,1],camera.R[2,2])
		#ry = atan2(-camera.R[2,0],sqrt(camera.R[2,1]**2 + camera.R[2,2]**2))
		#rz = atan2(camera.R[1,0], camera.R[0,0])
		#print(rx)
		#realR = realR*180/np.pi
		#actor.SetOrientation(rx, ry, rz)

		renderer.AddActor(actor)

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
			pod = Pod(folderPath + fileName,'',image_scale=1.0,to_meters=False,bg_image_filename='')
			cameras.append(pod.camera())

	window = MainWindow(cameras, vertices, faces)
	sys.exit(app.exec_())

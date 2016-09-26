import sys

import vtk
from PyQt5 import QtGui, QtCore, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from src.OBJIO import loadOBJ, writeOBJ

class MainWindow(QtWidgets.QMainWindow):
	def __init__(self, parent = None):
		QtWidgets.QMainWindow.__init__(self,parent)

		tabWidget = QtWidgets.QTabWidget();
		camVisibilityViewer = QVTKRenderWindowInteractor(tabWidget)
		renderer = vtk.vtkRenderer()
		camVisibilityViewer.GetRenderWindow().AddRenderer(renderer)
		iren = camVisibilityViewer.GetRenderWindow().GetInteractor()
		
		tabWidget.addTab(camVisibilityViewer, 'Camera Visibility')

		source = vtk.vtkSphereSource()
		source.SetCenter(0,0,0)
		source.SetRadius(5.0)

		mapper = vtk.vtkPolyDataMapper()
		mapper.SetInputConnection(source.GetOutputPort())

		actor = vtk.vtkActor()
		actor.SetMapper(mapper)
		renderer.AddActor(actor)

		self.setCentralWidget(tabWidget)
		self.show()
		iren.Initialize()

if __name__ == '__main__':
	app = QtWidgets.QApplication(sys.argv)

	vertices,faces,normals = loadOBJ(sys.argv[1])

	window = MainWindow()
	sys.exit(app.exec_())

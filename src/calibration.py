import os
import numpy as np
from src.thirdparty.body.loaders.scanner_frame import Pod
from opendr.camera import ProjectPoints

def cameraDistance(camPos1, camPos2, position):
	import numpy as np
	dist = np.linalg.norm(np.cross((position - camPos1),(position - camPos2))) / np.linalg.norm(camPos2 - camPos1)
	return dist

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
	def generateVisibilityMatrix(self, gridSize, gridScale):
		self.visibilityMatrix = np.zeros((gridSize[0], gridSize[1], gridSize[2]), dtype=np.uint8)

		for i in range(gridSize[0]):
			for j in range(gridSize[1]):
				for k in range(gridSize[2]):
					x = (i*gridScale[0]) - (gridSize[0]*gridScale[0]/2) + gridScale[0]/2
					y = (j*gridScale[1]) - (gridSize[1]*gridScale[1]/2) + gridScale[1]/2
					z = (k*gridScale[2]) - (gridSize[2]*gridScale[2]/2) + gridScale[2]/2

					self.ppointsA.v = [x,y,z]
					self.ppointsB.v = [x,y,z]

					x1 = self.ppointsA.r[0]
					y1 = self.ppointsA.r[1]
					x2 = self.ppointsB.r[0]
					y2 = self.ppointsB.r[1]

					if x1 < 1600 and x1 >= 0 and x2 < 1600 and x2 > 0 and y1 < 1200 and y1 >= 0 and y2 < 1200 and y2 >= 0:
						self.visibilityMatrix[k,j,i] = 1

def getStereoCamerasFromCalibration(folderPath):
	files = os.listdir(folderPath)
	cameras = []
	i = 0
	for fileName in files:
		if fileName.endswith('.tka'):
			pod = Pod(folderPath + fileName,'',image_scale=1.0,to_meters=False,bg_image_filename='')
			cameras.append((fileName[0:-4],pod.camera()))
			i += 1

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
	
	return stereoCameras

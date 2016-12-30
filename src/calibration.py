import os
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

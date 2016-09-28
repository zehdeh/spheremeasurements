#! /usr/bin/env python
__all__ = ["Camera", "Project"]

import numpy as np
import cv2
#import cv
import scipy.sparse as sp

class Camera(object):
    """Perspective projection camera model
    Follows the opencv distortion model.

    Attributes:
        camera_matrix: 3x3 camera matrix
        projection_matrix: 3x4 projection matrix
        f: 2x1 focal length
        c: 2x1 principle point
        k: 5x1 distortion params
        alpha: scalar skew coefficient
        R: 3x3 extrinsic rotation matrix
        r: 3x1 extrinsic rotation rodrigues vector
        T: 3x1 extrinsic translation
        w: scalar image width
        h: scalar image height
        sensor_size: (w,h) physical size of the sensor
        name: optional string with the camera name
    """
    def __init__(self, size=(640, 480), camera_matrix=None, k=None, r=None, R=None, t=None, f=None, c=None, alpha=None, sensor_size=(36,24), name=None, lookat=None):
        self.name = name
        self.w, self.h = size
        self.sensor_size = sensor_size
        
        if camera_matrix is None:
            if f is None:
                f = [500.,500.]
            if c is None:
                c = [self.w/2,self.h/2]
            if alpha == None:
                alpha = 0
            self.camera_matrix = np.array([[f[0], f[0]*alpha, c[0]],[0., f[1], c[1]],[0.,0.,1.]], dtype=np.float64)
        else:
            self.camera_matrix = np.require(camera_matrix, dtype=np.float64)
            
        if k is None:
            self.k = np.array([[0.,0.,0.,0.,0.]], dtype=np.float64).T
        else:
            self.k = np.require(k, dtype=np.float64)
        
        if t is None:
            self.t = np.array([[0.,0.,0.]], dtype=np.float64).T
        else:
            self.t = np.require(t, dtype=np.float64)

        if r is not None:
            self.r = np.require(r, dtype=np.float64)
        elif R is not None:
            self.r,J = cv2.Rodrigues(np.require(R, dtype=np.float64))
        else:
            self.r = np.array([[0.,0.,0.]], dtype=np.float64).T
    
        if lookat is not None :
            self.lookat(lookat)

    @property
    def R(self):
        R,J = cv2.Rodrigues(self.r)
        return R
    @property
    def f(self):
        return np.array([[self.camera_matrix[0,0], self.camera_matrix[1,1]]]).T
    @property
    def c(self):
        return np.array([[self.camera_matrix[0,2], self.camera_matrix[1,2]]]).T
    @property
    def alpha(self):
        return self.camera_matrix[0,1] / self.camera_matrix[0,0]
    @property
    def view_matrix(self):
        return np.hstack((self.R,self.t))
    @property
    def projection_matrix(self):
        return np.dot(self.camera_matrix, self.view_matrix)
    @property
    def fov(self):
        fovx, fovy, focalLength, principalPoint, aspectRatio = cv2.calibrationMatrixValues(self.camera_matrix, (self.w,self.h), self.sensor_size[0], self.sensor_size[1])
        return (fovx,fovy)
    @property
    def focalLength(self):
        fovx, fovy, focalLength, principalPoint, aspectRatio = cv2.calibrationMatrixValues(self.camera_matrix, (self.w,self.h), self.sensor_size[0], self.sensor_size[1])
        return focalLength
    @property
    def principalPoint(self):
        fovx, fovy, focalLength, principalPoint, aspectRatio = cv2.calibrationMatrixValues(self.camera_matrix, (self.w,self.h), self.sensor_size[0], self.sensor_size[1])
        return principalPoint
    @property
    def aspectRatio(self):
        fovx, fovy, focalLength, principalPoint, aspectRatio = cv2.calibrationMatrixValues(self.camera_matrix, (self.w,self.h), self.sensor_size[0], self.sensor_size[1])
        return aspectRatio

    @property
    def undistorted_camera(self):
        return Camera(camera_matrix=self.camera_matrix, r=self.r, t=self.t)

    def translate(self, delta_t):
        self.t += delta_t
    
    def rotate(self, delta_r_or_R):
        delta_R = cv2.Rodrigues(np.require(delta_r_or_R.flatten(), dtype=np.float64))[0] if (len(delta_r_or_R.flatten()) == 3) else delta_r_or_R
        self.r = cv2.Rodrigues(np.matrix(self.R)*np.matrix(delta_R))[0]
    
    @property
    def target(self):
        (origin, rays, xy_rect) = self.ray_for_image_points(np.array([self.c.flatten()]))
        return rays
    
    def lookat(self, target, up=None):
        look = target - self.position
        look = look / np.linalg.norm(look)
        if up is None:
            old_up = np.dot(self.R, np.array([[0,1,0]]).T)
            up = np.cross(look.T  , np.cross(self.target, old_up.T)).T
        up = up / np.linalg.norm(up)
        s = np.cross(look.T, up.T).T
        s = s / np.linalg.norm(s)
        u = np.cross(s.T, look.T).T
        self.r,J = cv2.Rodrigues(np.require( np.hstack((s,u,-look)) , dtype=np.float64))

    def project_points(self, xyz):
        if len(xyz.flatten()) == 0 :
            return np.array([])
        transpose = (xyz.shape[1] != 3)
        if transpose:
            xyz = xyz.T
        xy,J = cv2.projectPoints(xyz, self.r, self.t, self.camera_matrix, self.k)
        if transpose:
            xy = xy.T
        return xy.squeeze()
    
    def project_verts(self, verts, return_xyz=False):
        xy_chained_object = VertexProjection(verts,  self.camera_matrix, self.R, self.t, return_xyz)
#        xy,J = cv2.projectPoints(verts.r.reshape(-1,3), self.r, self.t, self.camera_matrix, 0.0*self.k)
#        check = np.sum(np.array(np.squeeze(xy) - xy_chained_object.r.reshape(-1,2))**2)
        return xy_chained_object
    
    # This method is useful for when a mesh has been converted from millimeters to meters
    # and the camera position should be adjusted accordingly
    def rescale_world_units(self, scale_factor=0.001):
        self.t = scale_factor*self.t
    
    @property
    def origin(self): # in camera space
        #return np.dot(self.R.T, np.zeros((3,1)).flatten() - self.t)
        return np.dot(self.R.T, np.zeros((3,1)) - self.t)
    @property
    def position(self): # in world space
        return - self.t
    
    # This is used to facilitate vertex visibility testing.
    # Basically, it defines a window that gives the CCD of the camera relative to its origin.
    @property
    def sensor_axis(self) :
        # Undistort takes into account the camera center, so we actually want the ray to the corner of the images
        #axis_pixels = np.array([self.c.flatten() + [0, self.w/2.0], self.c.flatten() + [self.h/2.0, 0], self.c.flatten()], dtype=float)
        axis_pixels = np.array([[self.w/2., self.h], [self.w, self.h/2.], [self.w/2., self.h/2.]], dtype=float)
        origin, rays, rectified_pixels = self.ray_for_image_points(axis_pixels)
        rays = [ray/np.linalg.norm(ray) for ray in [rays[1] - rays[2]*np.dot(rays[2],rays[1]), rays[0] - rays[2]*np.dot(rays[2],rays[0]), rays[2]]]
        rays = [(self.w/2.0)*rays[0], (self.h/2.0)*rays[1], np.mean(self.f)*rays[2]]
        return np.array(rays)
    
    def ray_for_image_points(self, xy):
        if len(xy.flatten()) == 0 :
            return (self.origin, np.array([]), np.array([]))
        transpose = (xy.shape[1] != 2)
        if transpose:
            xy = xy.T
        n = xy.shape[0]
        xy = np.array([[xy.T[0].flatten()],[xy.T[1].flatten()]]).T.copy() # copy is due to an obsolete bug in the old bindings
        
        #xy_rect = cv.CreateMat(xy.shape[0], 1, cv.CV_64FC2)
        # Annoyingly, UndistortPoints hasn't been ported to cv2
        xy_rect = cv2.undistortPoints(xy, self.camera_matrix.copy(), self.k.copy())
        
        xy_rect = np.asarray(xy_rect, dtype=float)
        # An annoying special case for when xy is just a single point
        xy_rect = xy_rect.squeeze() if xy_rect.shape[0] > 1 else np.array([[xy_rect[0, 0, 0]], [xy_rect[0, 0, 1]]]).T
        xy_rect = np.hstack((xy_rect, np.ones((n,1)))).T
        
        origin = self.origin
        rays = np.dot(self.R.T, xy_rect - self.t) - origin
        rays = rays / np.sum(rays**2, axis=0)**(0.5)
        
        if not transpose:
            origin = origin.T
            rays = rays.T
        return (origin, rays, xy_rect)
    
    def undistort_image(self, im):
        return cv2.undistort(im, self.camera_matrix, self.k)
    
    def redistort_image(self, im):
        raise Warning('Not implemented: redistort_image')
        return im

    @classmethod
    def calibrate(klass, xy, xyz, w, h):
        """Calibrate from one set of 2d-3d pairs. Not yet multi-image."""
        xy = np.require(xy, dtype='float32').T
        xyz = np.require(xyz, dtype='float32').T

        init_camera_matrix = np.zeros((3,3),'float32')
        init_camera_matrix[0,0]=500.
        init_camera_matrix[1,1]=500.
        init_camera_matrix[2,2]=1.
        init_camera_matrix[0,2]=w/2
        init_camera_matrix[1,2]=h/2
        init_dist_coefs = np.zeros(4,'float32')

        retval,camera_matrix,dist_coefs,rvecs,tvecs = cv2.calibrateCamera([xyz],[xy],(w,h),init_camera_matrix,init_dist_coefs, flags=cv2.CALIB_USE_INTRINSIC_GUESS)
        return klass((w,h), camera_matrix=camera_matrix, k=dist_coefs, r=rvecs[0], t=tvecs[0])

#!/usr/bin/env python
# encoding: utf-8
"""
scanner_frame.py

Created by Matthew Loper on 2013-02-28.
Copyright (c) 2013 MPI. All rights reserved.
"""

class Pod(object):
    def __init__(self, calib_fname, image_filename, image_scale, to_meters, bg_image_filename=None):
        from os.path import split, splitext
        self.calib_fname = calib_fname
        self.image_filename = image_filename
        self.bg_image_filename = bg_image_filename
        self.name = splitext(split(image_filename)[1])[0]
        self.id = self.name[:-1].zfill(2)
        self.to_meters = to_meters
        self.imcache = {}
        self.bgimcache = {}
        self.image_scale = image_scale

    def load_image(self, color=True, undistort=True):
        from cv2 import imread, CV_LOAD_IMAGE_COLOR, CV_LOAD_IMAGE_GRAYSCALE, resize
        from math import ceil
        from body.loaders.fcache import fc
        desired = (color, undistort)
        if desired not in self.imcache:
            im = imread(fc(self.image_filename), CV_LOAD_IMAGE_COLOR if color else CV_LOAD_IMAGE_GRAYSCALE)
            if self.bg_image_filename is not None:
                bg_im = imread(fc(self.bg_image_filename), CV_LOAD_IMAGE_COLOR if color else CV_LOAD_IMAGE_GRAYSCALE)

            if self.image_scale != 1.0:
                im = resize(im, (int(ceil(self.camera.w)), int(ceil(self.camera.h))))
                if self.bg_image_filename is not None:
                    bg_im = resize(bg_im, (int(ceil(self.camera.w)), int(ceil(self.camera.h))))
            if undistort:
                im = self.camera.undistort_image(im)
                if self.bg_image_filename is not None:
                    bg_im = self.camera.undistort_image(bg_im)
            self.imcache[desired] = im
            if self.bg_image_filename is not None:
                self.bgimcache[desired] = bg_im
                self.bg_im = self.bgimcache[desired]
        return self.imcache[desired]

    # Release cached stuff 
    def release(self):
        del(self.imcache)
        del(self.bgimcache)
        self.imcache = {}
        self.bgimcache = {}
        self._cached_properties.clear()

    #@cached_property
    def image(self):
        return self.load_image()


    #@cached_property
    def camera(self):
        # Turn the "%varname values" into a dict
        from re import findall, MULTILINE
        import numpy as np
        parms_raw = findall('%(\w+)\s*([\d\.\e\-\s]*)', open(self.calib_fname).read(), MULTILINE)
        parms = {p[0]: p[1].split() for p in parms_raw}
        parms = {pk: [float(v) for v in pv] for pk, pv in parms.items()}
        parms = {pk: pv[0] if len(pv)==1 else pv for pk, pv in parms.items()}
        c = parms

        R_pod = np.array(c['M']).reshape((3,3))
        c['is'] = np.asarray(np.array(c['is']) *  self.image_scale, np.int32)
        width = c['is'][0]
        height = c['is'][1]

        L_pod = np.array([c['X'], c['Y'], c['Z']])
        fl = c['f']
        px = c['x']
        py = c['y']
        cx = c['a']
        cy = c['b']
        k1 = c['K']
        k2 = c['K2']

        # dist parameters
        k1, k2 = self._convert_distortion_parms(k1, k2, fl, fl/px, fl/py, width, height)
        # k2 should be zero?
        dist = np.array([[k1[0], k2[0], 0, 0, 0]])

        fl *= self.image_scale
        cx *= self.image_scale
        cy *= self.image_scale

        fx = fl / px
        fy = fl / py

        camM = np.array([
            [fx,  0, cx ],
            [0,  fy, cy ],
            [0,   0, 1.0]]
        )

        #k1, k2 = self._convert_distortion_parms(k1, k2, fl_old, fx_old, fy_old, width, height)
        #dist = np.array([[k1[0], k2[0], 0, 0, 0]])

        T_pod = -R_pod.dot(L_pod)
        if self.to_meters:
            T_pod /= 1000.

        from os.path import splitext, basename
        from ..calib.camera import Camera
        from ..matlab.matlab import col
        name = splitext(basename(self.calib_fname))[0]
        return Camera(size=np.asarray(c['is'], np.int32), camera_matrix=camM, k=dist, R=R_pod, t=col(T_pod), name=name)


    def _convert_distortion_parms(self, k1, k2, fl, fx, fy, width, height):
        from scipy.linalg import lstsq
        import numpy as np
        # OpenCV wants radial distortion parameters that are applied to image plane coordinates
        # prior to being scaled by fx and fy (so not pixel coordinates). In contrast, k1 and k2
        # are defined via Tsai camera calibration http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/DIAS1/
        K1 = k1*(fl**2.0)
        K2 = k2*(fl**4.0)
        # Also, K1 and K2 are actually undistortion coefficients. They go from distorted to undistorted image
        # plane coordinates. OpenCV wants coefficients that go the other way.
        r_values = .01*np.array(range(1,101))*(((width/fx)**2.0 + (height/fy)**2.0)**0.5)
        undistorted_r_values = r_values*(1 + K1*(r_values**2.0) + K2*(r_values**4.0))
        distortion_factors = r_values/undistorted_r_values
        # Given the undistorted and distorted distances, we solve for the new distortion factors via linear regression
        k1, k2 = lstsq(np.matrix([undistorted_r_values**2.0, undistorted_r_values**4.0]).T, np.matrix(distortion_factors - 1.0).T)[0]
        return (k1, k2)


    def change_im_scale(self, new_image_scale):
        self.imcache = {}
        self.bgimcache = {}
        if hasattr(self, 'bg_img'):
            del bg_img
        if hasattr(self, '_cached_properties') and 'camera' in self._cached_properties:
            del self._cached_properties['camera']
        self.image_scale = new_image_scale



class ScannerFrame(object):
    def __init__(self, directory, image_scale=1.0, to_meters=True, scan_basename='*.obj'):
        if directory is not None:
            from os.path import join
            from body.loaders.fcache import remote_glob
            color_fnames = remote_glob(join(directory, '*C.png'))
            #calib_fnames = remote_glob(join(directory, 'calib_*C.tka'))
            from re import sub
            calib_fnames = [sub('/(\d+C)\.png$', '/calib_\\1.tka', p) for p in color_fnames]

            self.to_meters = True
            #self.scan_fname = remote_glob(join(directory, '*.obj'))[0]

            scan_fnames = remote_glob(join(directory, scan_basename))
            if scan_fnames:
                self.scan_fname = scan_fnames[0]
            else:
                self.scan_fname = remote_glob(join(directory, '*.obj'))[0]

            #self.pods = [Pod(calib_fnames[i], color_fnames[i], image_scale=image_scale, to_meters=to_meters) for i in range(len(calib_fnames))]
            self.pods = sorted([Pod(calib_fnames[i], color_fnames[i], image_scale=image_scale, to_meters=to_meters) for i in range(len(calib_fnames))], \
                key=lambda pod: pod.id)


    def set(self, scan_fname, pods, to_meters):
        self.scan_fname = scan_fname
        self.pods = pods
        self.to_meters = to_meters

    # Release cached images for memory
    def release(self):
        for pod in self.pods:
            pod.release()
        self._cached_properties.clear()    
 

    def mesh(self):
        from body.mesh import Mesh
        m = Mesh(filename = self.scan_fname)
        m.filename = self.scan_fname
        if self.to_meters:
            m.v /= 1000.
        return m


    def show(self, pod_indices=None):
        import numpy as np
        from body.mesh.meshviewer import MeshViewer
        from body.mesh.lines import Lines
        mv = MeshViewer()

        e = np.array([[0,1],[0,2],[0,3],[0,4]])
        all_lines = []
        if pod_indices == None:
            pod_indices = np.arange(len(self.pods))
        for idx in pod_indices:
            pod = self.pods[idx]
            print(pod.camera.fov)

            if True:
                import matplotlib.pyplot as plt
                plt.imshow(pod.load_image())
                plt.show()
            cam = pod.camera

            rays = cam.ray_for_image_points(np.array([
                0., 0.,
                cam.w, 0.,
                cam.w, cam.h,
                0., cam.h]).reshape((-1,2)))

            v = np.vstack((rays[0], rays[1]+rays[0]))

            all_lines.append(Lines(e=e, v=v))
        mv.set_static_lines(all_lines)
        mv.set_static_meshes([self.mesh])


    def compute_vertex_visibility(self, mesh=None):
        import numpy as np
        if mesh is None : mesh = self.mesh
        from body.mesh.visibility import visibility_compute
        if not hasattr(mesh, 'vn'):
            mesh.vn = mesh.estimate_vertex_normals()

        mesh.v = mesh.v.reshape(-1,3)
        vis,n_dot_cam = visibility_compute(v=mesh.v.reshape(-1,3), f=mesh.f, n=mesh.vn,
                                      cams=(np.array([ pod.camera.origin.flatten() for pod in self.pods ])),
                                      sensors=np.array([ axes.flatten() for axes in [ pod.camera.sensor_axis for pod in self.pods ]]))
        return (vis, n_dot_cam)


    def project_mesh_vertices(self, mesh=None, undistort=True):
        from cv2 import projectPoints
        from numpy import zeros
        if mesh is None : mesh = self.mesh
        projected_vertices = []
        for pod in self.pods:
            distCoeffs = zeros(5) if undistort else pod.camera.k
            (proj, J) = projectPoints(mesh.v, pod.camera.r, pod.camera.t, pod.camera.camera_matrix, distCoeffs=distCoeffs)
            projected_vertices.append(proj.squeeze())
        return projected_vertices


    def change_im_scale(self, new_image_scale):
        for pod in self.pods:
            pod.change_im_scale(new_image_scale)


def main():
    from body.images.offscreen_renderer import OffScreenRenderer as OSR
    from cv2 import imwrite, imshow
    for scale in [0.2, 0.4, 0.8, 1.0, 1.2, 2.0]:
        print("scale: %s" % scale)
        frame = ScannerFrame('/is/ps/shared/data/body/trials/models/50001_120921/50001_scans/CAESAR/', image_scale=scale)
        mesh = frame.mesh
        for i, pod in enumerate(frame.pods):
            im = pod.load_image(undistort=True)
            cam = pod.camera
            silh = OSR.render_silhouette(mesh, cam)
            im[:,:,2] = silh
            imwrite('/is/ps/shared/data/body/textured_alignments/scanner_frame_test/mine/sc'+str(scale)+'cam'+str(i)+'.png', im)
    exit()


    frame.show(pod_indices=[8])
    print(frame.pods[0].camera)
    imshow('abc', frame.pods[0].load_image())
    import pdb; pdb.set_trace()


if __name__ == '__main__':
    main()

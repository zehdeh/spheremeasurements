# spheremeasurements

* calibrations/ - contains Tsai calibration files used for different scans
* cleanup_mesh.py - used to isolate a sphere from a mesh with other objects in it
* dumbbell.py - a script to measure the distance between two fitted spheres. Not used yet.
* eliminate_badscans.py - Deletes files according to several criteria (0 vertices, extreme fitted radius, etc.)
* error-grid.py - The interactive visualization. You must provide a folder containing OBJ-files and a calibration folder
* focusmodel.py - Builds a 3d grid of expected errors by fitting a curve from the camera-focus experiment to scans
* generate_report.py - Generates a spreadsheet, or a csv file of measurement data. Needs a folder of OBJ-files
* import.sh - Some internal script used to copy data and clean it for all cameras
* mtf.py - Provided with a mesh with a 90deg edge, calculates the MTF
* plot_grid2d.py - Takes a 3D matrix and averages in one specified axis to view it as a 2d plot
* plot.py - Plots the gaussian- and mean curvature and the fitting error of a sphere mesh
* plot_spherefit.py - Plots points of a sphere mesh plus the fitted center point
* reconstructsh.py - Some attempt to calculate SH frequencies and reconstruct based on them. Not used yet.
* reports - Some older spreadsheets on experiments
* res/ - contains all the meshes
* src/ - commonly used code
* testsh.sh - calculates the sh frequencies for a given sphere mesh

Requirements:
-------------
Via apt-get:
* python-numpy
* python-scipy
* python-pyqt5
* python-matplotlib
* python-pip

Via pip:
* setuptools==12.0.5
* opendr

custom built:
* VTK7.1+

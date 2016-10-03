import sys
import numpy as np
from math import isnan
from scipy.special import sph_harm

# Create a sphere
r = 0.3
pi = np.pi
cos = np.cos
sin = np.sin
phi = np.linspace(0,pi,101)
theta = np.linspace(-pi,pi,101)

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

Lmax = 4
n = len(phi)
k = (Lmax + 1)**2
Y = np.zeros((n,k))



# Represent spherical harmonics on the surface of the sphere
for i in range(n):
	for l in range(Lmax+1):
		for m in range(-l,l+1):
			j = l**2 + l + m + 1
			s = sph_harm(l, m, theta[i], phi[i]).real
			#print 'l ' + str(l)
			#print 'm ' + str(m)
			if not isnan(s):
				Y[i,j-1] = s
print Y

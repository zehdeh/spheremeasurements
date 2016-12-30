import numpy as np
import cv2
from src.utils import rodrigues

M = np.array([[ 0.61967552, -0.0314627,  -0.78422723],[ 0.26199948, -0.93358761,  0.24447997],[-0.73983682, -0.35696538, -0.57027817]])
print 'M:'
print M

print 'Own result:'
print rodrigues(M)
print 'cv2 result:'
ans, J = cv2.Rodrigues(M)
print ans

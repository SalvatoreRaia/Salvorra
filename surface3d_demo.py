'''
======================
3D surface (color map)
======================

Demonstrates plotting a 3D surface colored with the coolwarm color map.
The surface is made opaque by using antialiased=False.

Also demonstrates using the LinearLocator and custom formatting for the
z axis tick labels.
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


figa, ax = plt.subplots(nrows=1, ncols=1, figsize=(10.80, 10.80))

x = np.arange(0, 10, 0.0000001)
y = (3*1000* 8.5**2* np.pi* (8.85*10**-12))

plt.show()
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
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
# coeffs
a = [0,1,1]
b = [2,2,2]
c = [3,3,0]
d = [10,10,10]

guess = [-7,5,5]
guessHist = [guess]
threshold = 0.1
alpha = 0.01

# error
errorHist = [];
def errorDef(d,guess):
    error = 0
    for i,num in enumerate(guess):
        error += (a[i]*checkIndex(guess,i-1) + b[i]*checkIndex(guess,i) + c[i]*checkIndex(guess,i+1) - d[i])**2
    return error

def checkIndex(guess,i):
    if i < 0 or i >= len(guess):
        return 0
    else:
        return guess[i]

X = np.arange(-7, -3, 0.25)
Y = np.arange(4, 8, 0.25)
X, Y = np.meshgrid(X, Y)
Z = 1.33
error = (X*b[0] + Y*c[0] - d[0])**2 + (X*a[1] + Y*b[1] + Z*c[1] - d[1])**2 + (Y*a[2] + Z*b[2] - d[2])**2

# Plot the surface.
surf = ax.plot_surface(X, Y, error, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0, 200)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

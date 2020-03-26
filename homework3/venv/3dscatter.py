import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

x = np.arange(-11, -1, 0.5)
y = np.arange(1,11,0.5)
z = np.arange(-5,5,0.5)

iterX = []
iterY = []
iterZ = []
for i in range(20):
    for j in range(20):
        for k in range(20):
            iterX.append(x[i])
            iterY.append(y[j])
            iterZ.append(z[k])

color = []
for i in range(len(iterX)):
    color.append((iterX[i] * 2 + iterY[i] * 3 - 10) ** 2 + (iterX[i] * 1 + iterY[i] * 2 + iterZ[i] * 3 - 10) ** 2 + (iterY[i] * 1 + iterZ[i] * 2 - 10) ** 2)

fig = plt.figure(figsize=(6, 6))
plt.scatter(iterX, iterY, iterZ,
           linewidths=1, alpha=.7,
           edgecolor='k',
           s = 40,

           c=color)
plt.colorbar()
plt.show()
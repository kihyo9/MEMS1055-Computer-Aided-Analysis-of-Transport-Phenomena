# Central differencing (2nd order) with Thomas method and gradient descent
# this works specifically for a tri-diagonal matrix, with error defined in hp.errorDef1
import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

#constants
k = 1
L = 10
sigma = 5.67e-8
h = 0.05
T_inf = 200
T_0 = 300
T_L = 400
C = (-1 / k) * (h * T_inf + sigma * T_inf ** 4)  # -100.72

# numerical solution parameters
nodes = 20  # grid points between BCs
dx = 10/(nodes+1)
initialGuess = [350]*nodes

# creation of coefficient functions and their partial derivatives
a = []
b = []
c = []
d = []

# __Change these functions to change the system__
# these are the coefficients - the terms are multiplied by x[i]
afuncCoeff = lambda x,i: 1
bfuncCoeff = lambda x,i: 2 + (h / k) * dx ** 2 + (sigma / k) * (dx ** 2) * (x[i] ** 3)
cfuncCoeff = lambda x,i: 1

dafunc = lambda x,i: 1
dbfunc = lambda x,i: 2 + (h / k) * dx ** 2 + 4*(sigma / k)*(dx**2)*x[i]**3
dcfunc = lambda x,i: 1

# __Coefficient Matrix and BCs__
for i in range(nodes):
    a.append(afuncCoeff)
    b.append(bfuncCoeff)
    c.append(cfuncCoeff)
    d.append(C*dx**2)
a[0] = lambda x,i: 0
c[-1] = lambda x,i: 0
d[0] = C*dx**2 - T_0
d[-1] = C*dx**2 - T_L

# iterate
errorHistory = []
iterHistory = []
coeffValues = hp.calculateCoeffs([a,b,c],initialGuess)
answer = hp.thomasAlgorithmSolver(coeffValues + [d])
iterHistory.append(answer)
errorHistory.append(hp.errorDef1(answer,initialGuess))

for i in range(100):
    coeffValues = hp.calculateCoeffs([a,b,c], answer)
    answer = hp.thomasAlgorithmSolver(coeffValues + [d])
    errorHistory.append(hp.errorDef1(answer, iterHistory[-1]))
    iterHistory.append(answer)

plt.figure()
plt.plot(errorHistory)
plt.yscale('log')
plt.title("Error history")

plt.figure()
plt.plot(answer)
plt.title("final answer")

plt.figure()
plt.yscale('log')
for i in range(len(list(zip(*iterHistory)))):
    plt.plot(list(zip(*iterHistory))[i])
plt.title("points")
plt.show()
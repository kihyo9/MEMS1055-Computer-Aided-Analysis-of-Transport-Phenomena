# Central differencing (2nd order) with Thomas method and gradient descent
# this works specifically for a tri-diagonal matrix, with error defined in hp.errorDef1
import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

#constants
k = 1
L = 10
sigma = 2.7e-9
h = 0.05
T_inf = 200
T_0 = 300
T_L = 400
C = (-1 / k) * (h * T_inf + sigma * T_inf ** 4)  # -14.32

# numerical solution parameters
nodes = 100  # grid points between BCs
dx = 10/(nodes+1)
initialGuess = [1]*nodes

# creation of coefficient functions and their partial derivatives
a = []
b = []
c = []
d = []

'''
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
'''

# __Change these functions to change the system__
# these are the coefficients - the terms are multiplied by x[i]
afuncCoeff = lambda x,i: 1
bfuncCoeff = lambda x,i: 2
cfuncCoeff = lambda x,i: 3

dafunc = lambda x,i: 1
dbfunc = lambda x,i: 2
dcfunc = lambda x,i: 3

# __Coefficient Matrix and BCs__
for i in range(nodes):
    a.append(afuncCoeff)
    b.append(bfuncCoeff)
    c.append(cfuncCoeff)
    d.append(10)
a[0] = lambda x,i: 0
c[-1] = lambda x,i: 0

# __partial derivative generation__
def checkIndex(val, index):
    if index < 0 or index >= nodes:
        return 0
    else:
        return val[i]


partials = []
f = lambda coeffs,x,i: coeffs[0][i](x,i)*checkIndex(x,i-1) + coeffs[1][i](x,i)*x[i] + coeffs[2][i](x,i)*checkIndex(x,i+1)
eq1 = lambda coeffs,x,i: coeffs[3][i] - f(coeffs,x,i)
deq1_dx1 = lambda coeffs,x,i: -1*dcfunc(x,i)
deq1_dx2 = lambda coeffs,x,i: -1*dbfunc(x,i)
deq1_dx3 = lambda coeffs,x,i: -1*dafunc(x,i)
for pNum in range(nodes):
    if pNum == 0:
        first_part2 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i))*deq1_dx2(coeffs,x,i)  # i
        first_part3 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i+1))*deq1_dx3(coeffs,x,i)  # i+1
        first_func = lambda coeffs,x,i: first_part2(coeffs,x,i)+first_part3(coeffs,x,i)
        partials.append(first_func)
    elif pNum == nodes - 1:
        last_part1 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i-1))*deq1_dx1(coeffs,x,i)  # i-1
        last_part2 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i))*deq1_dx2(coeffs,x,i)  # i
        last_func = lambda coeffs,x,i: last_part1(coeffs,x,i)+part2(coeffs,x,i)
        partials.append(last_func)
    else:
        part1 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i-1))*deq1_dx1(coeffs,x,i)  # i-1
        part2 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i))*deq1_dx2(coeffs,x,i)  # i
        part3 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i+1))*deq1_dx3(coeffs,x,i)  # i+1
        func = lambda coeffs,x,i: part1(coeffs,x,i)+part2(coeffs,x,i)+part3(coeffs,x,i)
        partials.append(func)


#####
alpha = 0.0001
errorThreshold = 0.0000001*nodes
answer = hp.gradientDescent2([a,b,c,d], initialGuess, hp.errorDef1, partials, alpha, errorThreshold)

plt.show()

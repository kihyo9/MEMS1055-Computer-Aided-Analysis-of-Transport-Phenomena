# Central differencing (2nd order) with Thomas method and gradient descent
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
nodes = 10  # grid points between BCs
# initialGuess = [0]*nodes
dx = 10/(nodes+1)
initialGuess = [350]*nodes

# initial guess -> first iteration of coeffs
a = []
b = []
c = []
d = []
for i in range(nodes):
    a.append(1)
    b.append(lambda x: 2 + (h/k)*dx**2 + (sigma/k)*(dx**2)*(x**3))
    c.append(1)
    d.append(C*dx**2)
a[0] = 0
c[-1] = 0
d[0] = C*dx**2 - T_0
d[-1] = C*dx**2 - T_L

# partial derivative generation
partials = []
# bfunc = lambda x: 2 + (h / k) * dx ** 2 + (sigma / k) * (dx ** 2) * (x ** 3)

def checkIndex(val, index):
    if index < 0 or index >= nodes:
        return 0
    else:
        return val[i]


f = lambda coeffs,x,i: coeffs[0][i]*checkIndex(x,i-1) + coeffs[1][i](x[i])*x[i] + coeffs[2][i]*checkIndex(x,i+1)
eq1 = lambda coeffs,x,i: x[i]- f(coeffs,x,i) - d[i]
for pNum in range(nodes):
    if pNum == 0:
        fist_part2 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i))*(1-3*(sigma / k)*(dx**2)*x[i]**2)  # i
        first_part3 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i+1))*(-1*coeffs[0][i])  # i+1
        first_func = lambda coeffs,x,i: fist_part2(coeffs,x,i)+first_part3(coeffs,x,i)
        partials.append(first_func)
    elif pNum == nodes - 1:
        last_part1 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i-1))*(-1*coeffs[2][i])  # i-1
        last_part2 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i))*(1-3*(sigma / k)*(dx**2)*x[i]**2)  # i
        last_func = lambda coeffs,x,i: last_part1(coeffs,x,i)+part2(coeffs,x,i)
        partials.append(last_func)
    else:
        part1 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i-1))*(-1*coeffs[2][i])  # i-1
        part2 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i))*(1-3*(sigma / k)*(dx**2)*x[i]**2)  # i
        part3 = lambda coeffs,x,i: 2*(eq1(coeffs, x, i+1))*(-1*coeffs[0][i])  # i+1
        func = lambda coeffs,x,i: part1(coeffs,x,i)+part2(coeffs,x,i)+part3(coeffs,x,i)
        partials.append(func)


#####
alpha = 0.001
errorThreshold = 0.0000001*nodes
hp.gradientDescent2([a,b,c,d], initialGuess, hp.errorDef1, partials, alpha, errorThreshold)

# thomas algorithm -> new result

# calculate error

# while
# evaluate error w threshold ?-> exit

# calculate partials(guess, new result) -> calculate new guess

# get new coeffs

# thomas algorithm -> new result

# calculate error
# end while




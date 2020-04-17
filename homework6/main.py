import numpy as np
import matplotlib.pyplot as plt
import explicit
import implicit
import upwind
import solver
import sys
################################################################
# Setup
################################################################

# Parameters
u = 1.
gammas = [0.,0.05, 1.,100.]
L = 1.
dx = 0.01

'''
# only for stability
r = gamma*dt/dx**2 # r < 1/2
G = 2*gamma/dt - u^2 # G>= 0

# based off the third condition
Re_x = 2.*dx/(u*dt) # Re_x <= 2/cfl
compare_Re_x = 2./cfl
'''


# More parameters
t_last = L/u
dx_steps = int(round(L/dx, 0))+1 # 101 steps
dx_data = [round(i*dx, 2) for i in range(dx_steps)]

# Initial condition
def initialCondition(x,L):
    isIterable = True
    try:
        iterator = iter(x)
        isIterable = True
    except TypeError:
        isIterable = False

    if isIterable:
        return [initialFunction(a) for a in x]
    else:
        return initialFunction(x)

def initialFunction(x):
    if x <= L/4. or x >= 3*L/4.:
        return 0
    elif x > L/4. and x <= L/2.:
        return x*4./L + -1
    else:
        return -x*4./L + 3

t_0 = [initialCondition(i*dx, L) for i in range(dx_steps)]

################################################################
# Solving
################################################################

# Choose solver
choice = input('Explicit (1), implicit (2), upwind (3): ')
while choice not in ['1','2','3']:
    print('Invalid input.')
    choice = input('Explicit (1), implicit (2), upwind (3): ')

# Choose gamma value
goodValue = False
gamma = 0
while not goodValue:
    try:
        gamma = float(input('Gamma value: '))
        goodValue = True
    except Exception:
        print("Not a valid value.")


# Solve
if choice is '1':
    fileName = 'explicit-data.txt'
    r = 0.25
    if gamma == 0:
        dt = 2.5e-4
    else:
        dt = r*dx**2/gamma
    dt_steps = int(np.ceil(t_last/dt))+1
    dt_data = [i*dt for i in range(dt_steps)]

    solver.printStability(dx, dt, u, gamma)
    s = input("Continue? ")
    if s == 'n':
        sys.exit()

    explicit.solve(dt, dx, gamma, u, L, t_0,fileName, dx_steps, dt_steps)
    solver.plotting(fileName, dx_data, dt_data, gamma)
elif choice is '2':
    fileName = 'implicit-data.txt'
    r = 0.25
    if gamma == 0:
        dt = 2.5e-4
    else:
        dt = r*dx**2/gamma
    dt_steps = int(np.ceil(t_last/dt))+1
    dt_data = [i*dt for i in range(dt_steps)]

    solver.printStability(dx, dt, u, gamma)
    s = input("Continue? ")
    if s == 'n':
        sys.exit()

    implicit.solve(dt, dx, gamma, u, L, t_0,fileName, dx_steps, dt_steps)
    solver.plotting(fileName, dx_data, dt_data, gamma)
elif choice is '3':
    fileName = 'upwind-data.txt'
    r = 0.25
    if gamma == 0:
        dt = 2.5e-4
    else:
        dt = r*dx**2/gamma
    dt_steps = int(np.ceil(t_last/dt))+1
    dt_data = [i*dt for i in range(dt_steps)]

    solver.printStability(dx, dt, u, gamma)
    s = input("Continue? ")
    if s == 'n':
        sys.exit()

    upwind.solve(dt, dx, gamma, u, L, t_0,fileName, dx_steps, dt_steps)
    solver.plotting(fileName, dx_data, dt_data, gamma)
else:
    print('How did you get here?')


'''
'import ...' gets a module object with can be used to call its methods and its classes
'from A import B' gets a class object as B. You cannot directly call class methods from A - call with A.B.method()
default class methods require a 'self' argument
class methods with decorators @classmethod and @staticmethod do not need 'self'
class methods have their class objects as its first argument
calling methods default class methods via a class object requires a 'self' input
instances of a class from a class object can call the class methods without a 'self' argument
'myclass = MyClass()' means instance = classObj() <- the 'constructor' __init__ is called
@classmethod can access a calling subclass but not a calling instance
'''
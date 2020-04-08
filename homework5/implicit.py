import matplotlib.pyplot as plt
import test
import numpy as np

# Simulation params
dt = 0.00001
dr = 0.05
t_inf = 10.
R = 1.
L = R/2.
T_0 = 400.
T_inf = 300.
dr_steps = int(round(R/dr, 0))
dt_steps = int(round(t_inf/dt, 0))

# Graphing params - water fluid, iron cylinder
k = 5. # of cylinder
h = 30.
cp = 4.18 # of water
rho = 1000. # of water
alpha = 23 # of cylinder #k_fluid/(cp*rho)
Bi = h*L/k
# t = dt*i
# Fo = alpha*t/L

### save F value to file ###
f = open("Tdata-implicit.txt", "w")

### Simulation ###
last_step = []
legend = []

# Data
T = []
t_0_step = [T_0 for i in range(dr_steps)]

# define coefficients
A = lambda x: (alpha / dr ** 2.) * (1. - (1. / (2. * x)))
B = -(2. * alpha / (dr ** 2.) + 1. / dt)
C = lambda x: (alpha / dr ** 2.) * (1. + (1. / (2. * x)))
D = -(1. / dt)

for n in range(dt_steps):
    print(str(n) + " out of " + str(dt_steps - 1.) + " timesteps")
    t_step = []

    # matrix coefficients
    Alist = []
    Blist = []
    Clist = []
    Dlist = []

    # initial condition of temp distribution
    if(n == 0):
        t_step = t_0_step.copy()
    else:
        # second-order central-differencing
        for m in range(dr_steps):
            # first gridpoint (actually second, but the i = 0 is undefined)
            if m == 0:
                Alist.append(0)
                Blist.append(1)
                Clist.append(-1)
                Dlist.append(0)
            # last gridpoint
            elif m == dr_steps - 1:
                Alist.append(-k/(k+h))
                Blist.append(1)
                Clist.append(0)
                Dlist.append(h*T_inf/(k+h))
            # other gridpoints
            else:
                Alist.append(A(m))
                Blist.append(B)
                Clist.append(C(m))
                Dlist.append(D*last_step[m])

        # evaluate the system of equations
        t_step = test.thomasAlgorithmSolver([Alist,Blist,Clist,Dlist])

    # copy use in the next timestep
    last_step = t_step.copy()

    # only add the data every so often
    if n < 100 and n % 15 == 0:
        T.append(last_step)
        legend.append(str(round(dt * n, 5)) + "s")
    elif n < 400 and n % 80 == 0:
        T.append(last_step)
        legend.append(str(round(dt * n, 5)) + "s")
    elif n < 2000 and n % 400 == 0:
        T.append(last_step)
        legend.append(str(round(dt * n, 5)) + "s")
    elif n % 500 == 0:
        T.append(last_step)
        legend.append(str(round(dt * n, 5)) + "s")

    # write to file
    for i,num in enumerate(t_step):
        if i == len(t_step) - 1:
            f.write(str(round(num,5)) + "\n")
        else:
            f.write(str(round(num,5)) + ", ")

    # it's pretty much steady state after this many timesteps
    if n >= 3000:
        break

f.close()

### end simulation ###

plt.figure()
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
selected_timesteps = [0, 1, 2, 3]
for i in range(len(T)):
    plt.plot([j*dr for j in range(dr_steps)], T[i])
plt.legend(legend)
plt.title('Plots of selected timesteps (implicit)')
plt.legend(legend,loc='center left',bbox_to_anchor=(1, 0.5))
plt.ylabel('Temperature [K]')
plt.xlabel('r [m]')
plt.show()
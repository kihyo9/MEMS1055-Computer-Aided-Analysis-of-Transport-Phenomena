import matplotlib.pyplot as plt
import numpy as np

# Simulation params
dt = 0.0001
dr = 0.1
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

# Stability check
G = alpha*dt/dr**2

# Data
T = []
# Boundary conditions
t_0_step = [T_0 for i in range(dr_steps)]

def firstGridpoint(T0, T1, T2, m):
    term1 = (1 / (m * dr)) * ((T1 - T0) / dr)
    term2 = (-3 * T0 + 4 * T1 - T2) / (dr**2)
    return alpha * (term1 + term2) * dt + T0

def lastGridpoint(T_m_1):
    return (h*T_inf + k*T_m_1)/(h+k)


def nextGridpoint(Tm1, T0, T1, m):
    term1 = (1 / (m * dr)) * (T1 - Tm1) / (2 * dr)
    term2 = (Tm1 - 2 * T0 + T1) / (dr ** 2)
    return alpha * (term1 + term2) * dt + T0


last_step = []
### Simulation ###
for n in range(dt_steps):
    print(str(n) + " out of " + str(dt_steps - 1) + " timesteps")
    t_step = []

    # initial condition of temp distribution
    if(n == 0):
        last_step = t_0_step.copy()
        t_step = t_0_step.copy()
    else:
        # second-order central-differencing
        for m in range(dr_steps):
            # first gridpoint (actually second, but the i = 0 is undefined)
            if m == 0:
                continue
            if m == 1:
                nextstep = firstGridpoint(last_step[m+0], last_step[m+1], last_step[m+2], m)
                t_step.append(nextstep) # the undefined gridpoint (i = 0) will be set equal to i = 1
                t_step.append(nextstep)
            # last gridpoint
            elif m == dr_steps - 1:
                nextstep = lastGridpoint(t_step[-1])
                t_step.append(nextstep)
            # other gridpoints
            else:
                nextstep = nextGridpoint(last_step[m - 1], last_step[m], last_step[m + 1],m)
                t_step.append(nextstep)
    if n % 15 == 0:
        T.append(t_step)
    last_step = t_step.copy()

    if n > 200:
        break

### end simulation ###

plt.figure()
# Shrink current axis by 20%
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
selected_timesteps = [0, 1, 2, 3]
legend = []
# for i in selected_timesteps:
for i in range(len(T)):
    plt.plot([j*dr for j in range(dr_steps)], T[i])
    legend.append(str(round(dt*i*15, 5)) + "s")
plt.title('Plots of selected timesteps')
plt.legend(legend,loc='center left',bbox_to_anchor=(1, 0.5))
plt.ylabel('Temperature [K]')
plt.xlabel('r [m]')
plt.show()
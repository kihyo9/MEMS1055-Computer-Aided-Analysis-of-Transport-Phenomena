import matplotlib.pyplot as plt
import numpy as np

dr = 0.005
R = 0.25
dr_steps = int(round(R/dr, 0)) + 1
legend = []
T = []

def explicitSolver(k = 5, h = 30,fileName = "Tdata-explicit.py"):
    # Simulation params
    L = R/2.
    T_0 = 400.
    T_inf = 300.

    # Graphing params - water fluid, iron cylinder
    cp = 4.18 # of water
    rho = 1000. # of water
    alpha = 23e-6 # of cylinder #k_fluid/(cp*rho)
    Bi = h*L/k
    # t = dt*i
    # Fo = alpha*t/L

    # Stability check
    G = 0.20 # set this
    dt = G*dr**2/alpha

    ### save F value to file ###
    f = open(fileName, "w")

    ### Simulation ###
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

    # Data
    t_0_step = [T_0 for i in range(dr_steps)]
    last_step = t_0_step.copy()
    n = -1

    while(last_step[0] > T_inf + 1):
        n += 1
        # print("timestep " + str(n) + ": " + str(round(dt*n,1)) + "s")
        t_step = []

        # initial condition of temp distribution
        if(n == 0):
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

        # only add the data every so often
        # if n < 100 and n % 15 == 0:
        #     T.append(t_step)
        #     legend.append(str(round(dt * n, 5)) + "s")
        # elif n < 400 and n % 80 == 0:
        #     T.append(t_step)
        #     legend.append(str(round(dt * n, 5)) + "s")
        # elif n < 2000 and n % 400 == 0:
        #     T.append(t_step)
        #     legend.append(str(round(dt * n, 5)) + "s")
        # elif n % 500 == 0:
        #     T.append(t_step)
        #     legend.append(str(round(dt * n, 5)) + "s")

        if n % 500 == 0:
            T.append(t_step)
            legend.append(str(round(dt * n, 0)) + "s")

        last_step = t_step.copy()

        # write to file
        for i,num in enumerate(t_step):
            if i == len(t_step) - 1:
                f.write(str(round(num,5)) + "\n")
            else:
                f.write(str(round(num,5)) + ", ")

    print("Stability criterion: " + str(round(G,3)))
    print("Steady-state time: " + str(n*dt))
    print("Timestep size: " + str(round(dt,2)) + "s")
    print("Fourier number: " + str(alpha*n*dt/R**2))
    f.close()

    ### end simulation ###

if __name__ == "__main__":
    explicitSolver(8,1)
    plt.figure()
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    for T_element in T:
        plt.plot([j*dr for j in range(dr_steps)], T_element)
    plt.legend(legend)
    plt.title('Plots of selected timesteps (explicit)')
    plt.legend(legend,loc='center left',bbox_to_anchor=(1, 0.5))
    plt.ylabel('Temperature [K]')
    plt.xlabel('r [m]')
    plt.show()
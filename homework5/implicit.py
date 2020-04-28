import matplotlib.pyplot as plt
import test
import numpy as np

dr = 0.005
R = 0.25
dr_steps = int(round(R/dr, 0)) + 1
legend = []
T = []

def implicitSolver(k = 5, h = 30,fileName = "Tdata-implicit.py"):
    # Simulation params
    L = R/2.
    T_0 = 400.
    T_inf = 300.

    # Graphing params - water fluid, iron cylinder
    # h
    alpha = 23e-6 # of cylinder #k_fluid/(cp*rho)
    # t = dt*i
    # Fo = alpha*t/L

    # Stability check
    G = 0.20 # set this
    dt = G*dr**2/alpha

    ### save F value to file ###
    f = open(fileName, "w")

    ### Simulation ###
    # Data
    t_0_step = [T_0 for i in range(dr_steps)]
    last_step = t_0_step.copy()
    n = -1

    # define coefficients
    A = lambda x: (alpha / dr ** 2.) * (1. - (1. / (2. * x)))
    B = -(2. * alpha / (dr ** 2.) + 1. / dt)
    C = lambda x: (alpha / dr ** 2.) * (1. + (1. / (2. * x)))
    D = -(1. / dt)

    while(last_step[0] > T_inf + 1):
        n += 1
        # print("timestep " + str(n) + ": " + str(round(dt*n,1)) + "s")
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

        # copy use in the next timestep
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
    implicitSolver()
    plt.figure()
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    for T_element in T:
        plt.plot([j*dr for j in range(dr_steps)], T_element)
    plt.legend(legend)
    plt.title('Plots of selected timesteps (implicit)')
    plt.legend(legend,loc='center left',bbox_to_anchor=(1, 0.5))
    plt.ylabel('Temperature [K]')
    plt.xlabel('r [m]')
    plt.show()
import solver
import time

def solve(dt, dx, gamma, u, t_0,fileName, dx_steps, dt_steps):
    solverType = "Upwind"

    ### clear file ###
    f = open(fileName, "w")
    f.close()


    # Data
    last_step = t_0.copy()
    linesWritten = 0
    starttime = time.time()
    for n in range(dt_steps):
        # show progress
        if n % 1000 == 0:
            print("timestep " + str(n) + ": " + str(round(dt*n,5)) + "s")
        t_step = []

        # initial condition of temp distribution
        if(n == 0):
            t_step = t_0.copy()
        else:
            # second-order central-differencing
            for m in range(dx_steps):
                # first gridpoint
                if m == 0:
                    nextstep = nextGridpoint(last_step[m - 1], last_step[m], last_step[0], dx, dt, u, gamma)
                    t_step.append(nextstep)
                # last gridpoint
                elif m == dx_steps - 1:
                    nextstep = nextGridpoint(last_step[m - 1], last_step[m], last_step[0], dx, dt, u, gamma)
                    t_step.append(nextstep)
                # other gridpoints
                else:
                    nextstep = nextGridpoint(last_step[m - 1], last_step[m], last_step[m + 1], dx, dt, u, gamma)
                    t_step.append(nextstep)

        # write to file
        newLine = solver.writeToFile(fileName, t_step, n, linesWritten, dt_steps-1)
        if newLine:
            linesWritten += 1

        # update last step
        last_step = t_step.copy()

    # performance
    solver.performanceTime(starttime, solverType, gamma, dt_steps)

    # result info
    print('\nThat took {} seconds'.format(time.time() - starttime))
    print("Timestep size: " + str(dt) + "s")
    print("Number of timesteps: " + str(dt_steps))


### Simulation ###
def firstGridpoint(T0, T1, T2, T3, dx, dt, u, gamma):
    term1 = gamma*(2*T0 - 5*T1 + 4*T2 - T3)/dx**3
    term2 = u*(-3*T0 + 4*T1 - T2)/(2*dx)
    return (term1-term2)*dt + T0


def lastGridpoint(T_first):
    return T_first


def nextGridpoint(Tm1, T0, T1, dx, dt, u, gamma):
    coeffTm1 = gamma*dt/(dx**2)
    coeffT0 = 1 + u*dt/dx- 2*gamma*dt/dx**2
    coeffT1 = gamma*dt/dx**2 - u*dt/dx
    return coeffTm1 * Tm1 + coeffT0 * T0 + coeffT1 * T1
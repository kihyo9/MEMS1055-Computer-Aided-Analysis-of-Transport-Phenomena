import solver
import time
import multiprocessing
import poolTest


'''
Fix:
- the results of the pool processes are not being stored anywhere; when each process starts, the only variables available are the input parameters
- thus, the only iteration showing up is the initial condition
'''

### Simulation ###
def nextGridpoint(m, last_step, _dt, _dx, _gamma, _u, _dx_steps,t_step):
    # m, last_step, _dt, _dx, _gamma, _u, _dx_steps = data[0], data[1], data[2], data[3], data[4], data[5], data[6]
    try:
        # print(id(last_step))
        coeffTm1 = _u * _dt / (2 * _dx) + _gamma * _dt / (_dx ** 2)
        coeffT0 = 1 - 2 * _gamma * _dt / _dx ** 2
        coeffT1 = _gamma * _dt / _dx ** 2 - _u * _dt / (2 * _dx)
        t_step[m] = coeffTm1 * last_step[(_dx_steps + m - 1) % _dx_steps] + coeffT0 * last_step[m] + coeffT1 * last_step[(_dx_steps + m + 1) % _dx_steps]
        # print("Iteration " + str(m))
    except (IndexError) as e:
        print("Iteration " + str(m) + ", dx is " + str(_dx))
        print(last_step)
        print(last_step[(_dx_steps + m - 1) % _dx_steps])
        print(last_step[m])
        print(last_step[(_dx_steps + m + 1) % _dx_steps])
        print(_dx_steps)

def setGlobals(dt, dx, gamma, u,dx_steps):
    global _dx_steps
    _dx_steps = dx_steps
    global _gamma
    _gamma= gamma
    global _dt
    _dt= dt
    global _dx
    _dx= dx
    global _u
    _u= u

def solve(dt, dx, gamma, u, L, t_0,fileName, dx_steps, dt_steps):
    setGlobals(dt, dx, gamma, u,dx_steps)

    solverType = "Explicit2"

    ### clear file ###
    f = open(fileName, "w")
    f.close()

    # Data
    global last_step
    last_step = t_0.copy()
    global t_step
    t_step = [0 for _ in range(_dx_steps)]
    linesWritten = 0
    starttime = time.time()
    for n in range(dt_steps):
        # show progress
        if n % 1 == 0:
            print("timestep " + str(n) + ": " + str(round(dt*n,5)) + "s")

        # initial condition of temp distribution
        if(n == 0):
            t_step = t_0.copy()
        else:
            # second-order central-differencing
            pool = multiprocessing.Pool()
            pool.starmap(nextGridpoint, [(i,last_step, dt, dx, gamma, u,dx_steps, t_step) for i in range(dx_steps)])
            pool.close()


        # write to file
        newLine = solver.writeToFile(fileName, t_step, n, linesWritten, dt_steps-1)
        if newLine:
            linesWritten += 1

        # update last step
        last_step = t_step.copy()

    # performance
    solver.performanceTime(starttime, solverType, gamma, dt_steps,u,L)

    # result info
    print('\nThat took {} seconds'.format(time.time() - starttime))
    print("Timestep size: " + str(dt) + "s")
    print("Number of timesteps: " + str(dt_steps))


def basic_func(x):
    if x == 0:
        return 'zero'
    elif x%2 == 0:
        return 'even'
    else:
        return 'odd'

def multiprocessing_func(x):
    y = x*x
    time.sleep(2)
    print('{} squared results in a/an {} number'.format(x, basic_func(y)))

if __name__ == "__main__":
    # processes

    # processes = []
    # for i in range(0, 10):
    #     p = multiprocessing.Process(target=multiprocessing_func, args=(i,))
    #     processes.append(p)
    #     p.start()
    #
    # for process in processes:
    #     process.join()

    # pools
    poolTest.poolTest()
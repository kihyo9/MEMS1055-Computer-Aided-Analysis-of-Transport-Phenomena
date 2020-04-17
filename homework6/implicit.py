import solver
import time
import numpy as np

def solve(dt, dx, gamma, u, t_0,fileName, dx_steps, dt_steps):
    solverType = "Implicit"

    ### clear file ###
    f = open(fileName, "w")
    f.close()

    # define constant part of coefficients
    A = (gamma / dx ** 2.) + u/(2*dx)
    B = -2. * gamma / (dx ** 2.) + -1. / dt
    C = (gamma / dx ** 2.) - u / (2. * dx)
    D = -(1. / dt)

    # corners of the matrix
    # first and last gridpoints are adjacent
    alpha = C # bottom-left
    beta = A # top-right
    gamma2 = -B

    # Data
    last_step = t_0.copy()
    linesWritten = 0
    starttime = time.time()
    for n in range(dt_steps):

        # show progress
        if n % 1000 == 0:
            print("timestep " + str(n) + ": " + str(round(dt*n,5)) + "s")

        # matrix coefficients
        Alist = [A for _ in range(dx_steps)]
        Blist = [B for _ in range(dx_steps)]
        Clist = [C for _ in range(dx_steps)]
        Dlist = [D*last_step[i] for i in range(dx_steps)]

        if n == 0:
            # initial condition of temp distribution
            t_step = t_0.copy()
        else:
            # Modified Thomas algorithm
            t_step = thomasAlgorithmSolver2([Alist,Blist,Clist,Dlist],alpha, beta, gamma2)

        # write to file
        newLine = solver.writeToFile(fileName, t_step, n, linesWritten, dt_steps-1)
        if newLine:
            linesWritten += 1

        # current step is now previous timestep
        last_step = t_step.copy()

    # performance
    solver.performanceTime(starttime, solverType, gamma, dt_steps)

    # result info
    print('\nThat took {} seconds'.format(time.time() - starttime))
    print("Timestep size: " + str(dt) + "s")
    print("Number of timesteps: " + str(dt_steps))

def thomasAlgorithmSolver(matrix):
    a = matrix[0]
    b = matrix[1]
    c = matrix[2]
    d = matrix[3]

    n = len(d)

    cc = []
    dd = []
    for i in range(n):
        if i == 0:
            cc.append(c[i]/b[i])
            dd.append(d[i] / b[i])
        elif i == n - 1:
            dd.append((d[i] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1]))
        else:
            cc.append(c[i]/(b[i] - a[i]*cc[i-1]))
            dd.append((d[i] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1]))

    x = [0 for i in range(n)]
    for i in range(n):
        q = n-1-i

        if q == n-1:
            x[q] = dd[q]
        else:
            x[q] = dd[q] - cc[q]*x[q+1]
    return x

def thomasAlgorithmSolver2(matrix,alpha, beta, gamma2):
    a = matrix[0].copy()
    b = matrix[1].copy()
    c = matrix[2].copy()
    d = matrix[3].copy()

    n = len(d)

    # corrections
    b[0] -= gamma2
    b[-1] -= alpha*beta/gamma2

    # intermediate vectors
    u = [0. for _ in range(n)]
    u[0] = gamma2
    u[-1] = alpha

    v = [0. for _ in range(n)]
    v[0] = 1.
    v[-1] = beta/gamma2

    y = thomasAlgorithmSolver([a,b,c,d])
    z = thomasAlgorithmSolver([a,b,c,u])

    return vectorMath(v,y,z)

def vectorMath(v,y,z):
    term1 = sum([f*g for f,g in zip(v,y)]) # scalar
    term2 = 1. + sum([f*g for f,g in zip(v,z)]) # scalar
    return [f - (term1/term2)*g for f,g in zip(y,z)]


if __name__ == "__main__":
    n = 6
    a = [1. for _ in range(n)]
    b = [2. for _ in range(n)]
    c = [3. for _ in range(n)]
    d = [10. for _ in range(n)]


    alpha = 2.
    beta = 2.
    gamma = -1*b[0]

    a[0] = beta
    c[-1] = alpha

    print("Direct inverse")
    AA = np.array([[2,3,0,0,0,2],[1,2,3,0,0,0],[0,1,2,3,0,0],[0,0,1,2,3,0],[0,0,0,1,2,3],[2,0,0,0,1,2]])
    BB = np.array([10,10,10,10,10,10])
    direct = np.linalg.solve(AA,BB)
    print(direct)


    print("\nSherman Morrison")

    y = thomasAlgorithmSolver2([a,b,c,d],alpha, beta, gamma)
    for i,x in enumerate(y):
        if i == 0:
            print("%f, diff is %f" % (x, a[i]*y[-1]+b[i]*y[i]+c[i]*y[1] - d[i]))
        elif i == n - 1:
            print("%f, diff is %f" % (x, a[i]*y[i-1]+b[i]*y[i]+c[i]*y[0] - d[i]))
        else:
            print("%f, diff is %f" % (x, a[i]*y[i-1]+b[i]*y[i]+c[i]*y[i+1] - d[i]))

    print(np.allclose(np.dot(AA, y), BB))
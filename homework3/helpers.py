import numpy as np
import matplotlib.pyplot as plt

def analyticalSol(x, theta_start):
    y = []
    for i, num in enumerate(x):
        y.append(theta_start*np.cos(4.04351*x[i]))
    return y

def gradientDescent2(coeffs, initialGuess, errorDef, partials, alpha, errorThreshold):

    errorHistory = []
    iterHistory = []
    iterHistory.append(initialGuess)

    # calculate coeffs
    calculatedCoeffs = coeffs.copy()
    newB = []
    for i,func in enumerate(coeffs[1]):
        newB.append(func(initialGuess[i]))
    calculatedCoeffs[1] = newB

    # calculate new result
    newResult = triDiagSubSolver(calculatedCoeffs,initialGuess)

    # initialError
    error = errorDef(initialGuess, newResult)
    errorHistory.append(error)

    #first try!
    if(error < errorThreshold):
        return initialGuess
    count = 1
    print("Iteration: " + str(count) + ", error: " + str(error))

    #first new guess with partials
    newIter = []
    for i, partial in enumerate(partials):
        newIter.append(initialGuess[i] - alpha * partial(coeffs, initialGuess, i))
    iterHistory.append(newIter)

    #calculate coeffs
    newB = []
    for i,func in enumerate(coeffs[1]):
        newB.append(func(newIter[i]))
    calculatedCoeffs[1] = newB

    # calculate new result
    newResult = triDiagSubSolver(calculatedCoeffs,newIter)

    #calculate error
    error = errorDef(newIter, newResult)
    errorHistory.append(error)

    #iterations of gradient descent
    # while(error > errorThreshold):
    for i in range(250):
        count += 1
        print("Iteration: " + str(count) + ", error: " + str(error))
        oldIter = newIter.copy()

        # new guess with partials
        newIter = []
        for i, partial in enumerate(partials):
            newIter.append(oldIter[i] - alpha*partial(coeffs, oldIter, i))
        iterHistory.append(newIter)

        # newResult
        newB = []
        for i, func in enumerate(coeffs[1]):
            newB.append(func(newIter[i]))
        calculatedCoeffs[1] = newB
        newResult = triDiagSubSolver(calculatedCoeffs,newIter)

        # calculate new error
        error = errorDef(newIter, newResult)
        errorHistory.append(error)

    plt.figure()
    plt.title('Iteration history of grid points')
    y = zip(*iterHistory)
    z = list(y)
    for i in range(10):
        plt.plot(list(range(252)),z[i])
    plt.xlabel('iteration')
    plt.ylabel('Temperature [K]')

    plt.figure()
    plt.title('Error history of grid points')
    plt.plot(list(range(252)),errorHistory)
    plt.xlabel('iteration')
    plt.ylabel('Error')
    plt.show()


    return newIter

def errorDef1(initials, newIter):
    errorSum = 0
    for i, initial in enumerate(initials):
        errorSum = errorSum + (initials[i] - newIter[i])**2
    return errorSum

def thomasAlgorithmSolver(matrix):
    a = matrix[0]
    b = matrix[1]
    c = matrix[2]
    d = matrix[3]

    n = len(d)
    cc = []
    dd = []
    x = []

    # cc
    for i in range(n-1):
        if i == 1:
            cc.append(c[i]/b[i])
        else:
            cc.append(c[i] / (b[i] - a[i]*cc[i-1]))

    # dd
    for i in range(n):
        if i == 1:
            dd.append(d[i]/b[i])
        else:
            dd.append((d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]))

    # x
    for i in range(n):
        i = n - i

        if i == n:
            x.append(dd[n])
        else:
            x.append(dd[i]/(cc[i]*x[i+1]))

    return x

def soeSolver(input, soe):
    newResult = []
    for equation in soe:
        newResult.append(equation(input))
    return newResult

def triDiagSubSolver(coeffs, guess):
    x = []
    for i in range(len(guess)):
        if i == 0:
            x.append(coeffs[1][i]+coeffs[2][i])
        elif i == len(guess) - 1:
            x.append(coeffs[0][i]+coeffs[1][i])
        else:
            x.append(coeffs[0][i]+coeffs[1][i]+coeffs[2][i])

    return x

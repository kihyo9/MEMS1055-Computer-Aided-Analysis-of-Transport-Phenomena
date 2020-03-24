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
    calculatedCoeffs = [0,0,0]
    newA = []
    newB = []
    newC = []
    for i,func in enumerate(list(zip(*coeffs))):
        newA.append(func[0](initialGuess, i))
        newB.append(func[1](initialGuess, i))
        newC.append(func[2](initialGuess, i))
    calculatedCoeffs[0] = newA
    calculatedCoeffs[1] = newB
    calculatedCoeffs[2] = newC

    # calculate new result
    newResult = triDiagSubSolver(calculatedCoeffs,initialGuess)

    # initialError
    error = errorDef(coeffs[3], newResult)
    errorHistory.append(error)
    lowestError = [error,0]

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
    newA = []
    newB = []
    newC = []
    for i,func in enumerate(list(zip(*coeffs))):
        newA.append(func[0](newIter, i))
        newB.append(func[1](newIter, i))
        newC.append(func[2](newIter, i))
    calculatedCoeffs[0] = newA
    calculatedCoeffs[1] = newB
    calculatedCoeffs[2] = newC

    # calculate new result
    newResult = triDiagSubSolver(calculatedCoeffs,newIter)

    # calculate error
    error = errorDef(coeffs[3], newResult)
    errorHistory.append(error)
    if error < lowestError[0]:
        lowestError = [error,count]


    maxIterations = 1000

    # iterations of gradient descent
    # while(error > errorThreshold):
    for i in range(maxIterations):
        count += 1
        print("Iteration: " + str(count) + ", error: " + str(error))
        oldIter = newIter.copy()

        # new guess with partials
        newIter = []
        for i, partial in enumerate(partials):
            newIter.append(oldIter[i] - alpha*partial(coeffs, oldIter, i))
        iterHistory.append(newIter)

        # new coeffs
        calculatedCoeffs = coeffs.copy()
        newA = []
        newB = []
        newC = []
        for i, func in enumerate(list(zip(*coeffs))):
            newA.append(func[0](newIter, i))
            newB.append(func[1](newIter, i))
            newC.append(func[2](newIter, i))
        calculatedCoeffs[0] = newA
        calculatedCoeffs[1] = newB
        calculatedCoeffs[2] = newC

        # new result
        newResult = triDiagSubSolver(calculatedCoeffs,newIter)

        # calculate new error
        error = errorDef(coeffs[3], newResult)
        errorHistory.append(error)
        if error < lowestError[0]:
            lowestError = [error, count]

    plt.figure()
    plt.title('Iteration history of grid points')
    y = zip(*iterHistory)
    z = list(y)
    for i in range(len(z)):
        plt.plot(list(range(maxIterations+2)),z[i])
    plt.xlabel('iteration')
    plt.ylabel('Temperature [K]')

    plt.figure()
    plt.title('Error history of grid points')
    plt.plot(list(range(maxIterations+2)),errorHistory)
    plt.xlabel('iteration')
    plt.ylabel('Error')

    plt.figure()
    plt.title('Lowest error plot')
    plt.plot(np.linspace(0,10,len(iterHistory[lowestError[1]])),iterHistory[lowestError[1]])
    plt.xlabel('position [m]')
    plt.ylabel('Temperature [K]')
    print(lowestError[1])
    print(iterHistory[lowestError[1]])
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
            x.append(coeffs[1][i]*guess[i]+coeffs[2][i]*guess[i+1])
        elif i == len(guess) - 1:
            x.append(coeffs[0][i]*guess[i-1]+coeffs[1][i]*guess[i])
        else:
            x.append(coeffs[0][i]*guess[i-1]+coeffs[1][i]+coeffs[2][i]*guess[i+1])

    return x

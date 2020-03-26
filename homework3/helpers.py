import numpy as np
import matplotlib.pyplot as plt

def analyticalSol(x, theta_start):
    y = []
    for i, num in enumerate(x):
        y.append(theta_start*np.cos(4.04351*x[i]))
    return y

def gradientDescent2(coeffs, initialGuess, partials, alpha, errorThreshold):
    # iterations of guesses
    iterHistory = []
    shortiterHistory = []
    iterHistory.append(initialGuess)
    shortiterHistory.append(initialGuess)

    # history of errors
    errorHistory = []
    shortErrorHistory = []
    error = newGuessToError(coeffs,initialGuess)
    errorHistory.append(error)
    shortErrorHistory.append(error)

    #newIter variable
    newIter = initialGuess.copy()
    count = 0

    #while (error > errorThreshold):
    for num in range(132):
        # determine grid point with largest dError/dx
        largestPartial = [0, 0]
        for i, partial in enumerate(partials):
            currentPartial = partial(coeffs, initialGuess, i)
            if i == 0:
                largestPartial[0] = currentPartial
                largestPartial[1] = i
            else:
                if abs(largestPartial[0]) > abs(currentPartial):
                    largestPartial[0] = currentPartial
                    largestPartial[1] = i

        # find error of new guess
        newIter[largestPartial[1]] -= alpha * largestPartial[0]
        iterHistory.append(newIter)
        subError = newGuessToError(coeffs, newIter)
        errorHistory.append(subError)
        prevSubError = errorHistory[-1]

        while(1 - subError/prevSubError > 0.001):
            prevSubError = subError
            newIter[largestPartial[1]] -= alpha * largestPartial[0]
            iterHistory.append(newIter)
            subError = newGuessToError(coeffs, newIter)
            errorHistory.append(subError)

        count += 1
        print("Iteration " + str(count) + "... " + "Moving grid point: " + str(largestPartial[1]) + "... " + "Error: " + str(errorHistory[-1]))
        shortiterHistory.append(iterHistory[-1])
        shortErrorHistory.append(errorHistory[-1])

    x = []
    for i, partial in enumerate(partials):
        x.append(partial(coeffs, iterHistory[-1], i))

    plt.figure()
    plt.title('Iteration history of grid points')
    y = zip(*iterHistory)
    z = list(y)
    for i in range(len(z)):
        plt.plot(list(range(len(iterHistory))),z[i])
    plt.xlabel('iteration')
    plt.ylabel('Temperature [K]')

    plt.figure()
    plt.title('Error history of grid points')
    plt.plot(list(range(len(errorHistory))),errorHistory)
    plt.xlabel('iteration')
    plt.ylabel('Error')

    return newIter

def errorDef1(initials, newIter):
    errorSum = 0
    for i, initial in enumerate(initials):
        errorSum = errorSum + (initials[i] - newIter[i])**2
    return errorSum

def newGuessToError(coeffs, guess):
    # calculate coeffs
    calculatedCoeffs = calculateCoeffs(coeffs, guess)

    # calculate new result
    newResult = triDiagSubSolver(calculatedCoeffs,guess)

    # error
    error = errorDef1(coeffs[3], newResult)

    return error

def calculateCoeffs(coeffs, guess):
    newA = []
    newB = []
    newC = []
    for i,func in enumerate(list(zip(*coeffs))):
        newA.append(func[0](guess, i))
        newB.append(func[1](guess, i))
        newC.append(func[2](guess, i))
    return [newA, newB, newC]

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
        if i == 0:
            cc.append(c[i]/b[i])
        else:
            cc.append(c[i] / (b[i] - a[i]*cc[i-1]))

    # dd
    for i in range(n):
        if i == 0:
            dd.append(d[i]/b[i])
        else:
            dd.append((d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]))

    # x
    for i in range(n):
        i = n - i - 1

        if i == n-1:
            x.append(dd[i])
        else:
            x = [dd[i]/(cc[i]*x[0])] + x

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

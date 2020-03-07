import numpy as np

def analyticalSol(x, theta_start):
    y = []
    for i, num in enumerate(x):
        y.append(theta_start*np.cos(4.04351*x[i]))
    return y

def gradientDescent(soe, initialGuess, errorDef, partials, alpha, errorThreshold):

    #make initial newIter
    newResult = []
    for equation in soe:
        newResult.append(equation(initialGuess))

    # initialError
    error = errorDef(initialGuess, newResult)

    #first try!
    if(error < errorThreshold):
        return initialGuess

    #first new guess
    newIter = []
    for i, partial in enumerate(partials):
        newIter.append(initialGuess[i] - alpha * partial(initialGuess))

    newResult = []
    for equation in soe:
        newResult.append(equation(newIter))
    error = errorDef(newIter, newResult)

    #iterations of gradient descent
    while(error > errorThreshold):
        oldIter = newIter.copy()
        newIter = []

        for i, partial in enumerate(partials):
            newIter.append(oldIter[i] - alpha*partial(oldIter))

        newResult = []
        for equation in soe:
            newResult.append(equation(newIter))

        error = errorDef(newIter, newResult)

    return newIter

def errorDef1(initials, newIter):
    errorSum = 0
    for i, initial in enumerate(initials):
        errorSum = errorSum + (initials[i] - newIter[i])**2
    return errorSum
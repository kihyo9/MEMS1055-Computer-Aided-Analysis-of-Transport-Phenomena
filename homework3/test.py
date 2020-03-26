# temperature dist of the bar
# shooting+rk4 is finished
# finite difference with thomas method
# determine optimum dx for each case
## maybe find the difference in solution with each successive decrease in dx?

import numpy as np
import matplotlib.pyplot as plt

# coeffs
a = [0,1,1]
b = [2,2,2]
c = [3,3,0]
d = [10,10,10]

guess = [-7,5,5]
guessHist = [guess]
threshold = 0.1
alpha = 0.01

# error
errorHist = [];
def errorDef(d,guess):
    error = 0
    for i,num in enumerate(guess):
        error += (a[i]*checkIndex(guess,i-1) + b[i]*checkIndex(guess,i) + c[i]*checkIndex(guess,i+1) - d[i])**2
    return error

def checkIndex(guess,i):
    if i < 0 or i >= len(guess):
        return 0
    else:
        return guess[i]

error = errorDef(d,guess)
errorHist.append(error)
count = 0
for i in range(1000):
    row = count % 3
    # partial of dError/dx[1]
    if row == 0:
        partial = 2*(b[0]*checkIndex(guess,0) + c[0]*checkIndex(guess,1) - d[0])*b[0]
        partial += 2*(a[1]*checkIndex(guess,0) + b[1]*checkIndex(guess,1) + c[1]*checkIndex(guess,2) - d[1])*a[1]
    elif row == 1:
        partial = 2*(b[0]*checkIndex(guess,0) + c[0]*checkIndex(guess,1) - d[0])*c[0]
        partial += 2*(a[1]*checkIndex(guess,0) + b[1]*checkIndex(guess,1) + c[1]*checkIndex(guess,2) - d[1])*b[1]
        partial += 2 * (a[2] * checkIndex(guess, 1) + b[2] * checkIndex(guess, 2) + - d[2]) * a[2]
    elif row == 2:
        partial += 2*(a[1]*checkIndex(guess,0) + b[1]*checkIndex(guess,1) + c[1]*checkIndex(guess,2) - d[1])*b[2]
        partial += 2 * (a[2] * checkIndex(guess, 1) + b[2] * checkIndex(guess, 2) + - d[2]) * c[1]

    guess[row] -= alpha*partial
    copyGuess = guess.copy()
    guessHist.append(copyGuess)

    errorHist.append(errorDef(d,guess))

    count += 1

plt.figure()
plt.title('error plot')
plt.plot(errorHist)
plt.figure()
plt.title('guess[i] over iterations')
plt.plot(list(zip(*guessHist))[0])
plt.plot(list(zip(*guessHist))[1])
plt.plot(list(zip(*guessHist))[2])
plt.legend(['x[0]','x[1]','x[2]'])
plt.show()
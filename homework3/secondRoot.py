# second root

import matplotlib.pyplot as plt
import numpy as np

# rk4

k = 1
L = 10
sigma = 2.7e-9
h = 0.05
T_inf = 200
T_0 = 300
T_L = 400
C = (-1 / k) * (h * T_inf + sigma * T_inf ** 4)  # -14.32


# define rk4 calculation
def rk4(current_T, current_Tprime):
    T1 = current_Tprime * dx
    Tprime1 = Tprimeprime(current_T) * dx

    T2 = (current_Tprime + 0.5 * Tprime1) * dx
    Tprime2 = Tprimeprime(current_T + 0.5 * T1) * dx

    T3 = (current_Tprime + 0.5 * Tprime2) * dx
    Tprime3 = Tprimeprime(current_T + 0.5 * T2) * dx

    T4 = (current_Tprime + Tprime3) * dx
    Tprime4 = Tprimeprime(current_T + T3) * dx

    new_T = current_T + (1 / 6) * (T1 + 2 * T2 + 2 * T3 + T4)
    new_Tprime = current_Tprime + (1 / 6) * (Tprime1 + 2 * Tprime2 + 2 * Tprime3 + Tprime4)

    return new_T, new_Tprime


# define d2T/dx2
def Tprimeprime(current_T):
    return (h / k) * current_T + (sigma / k) * current_T ** 4 + C


# numerical solution parameters
dx = 0.001
leftEnd, rightEnd = -45, -40
sList = np.linspace(leftEnd, rightEnd, (rightEnd - leftEnd) * 4 + 1)  # Guesses for dt/dx @ x=0
T_LList = []  # List of errors from the exact BC: T(L) = 400

for s in sList:
    # initial conditions
    T = [T_0]
    Tprime = [s]  # dT/dx

    # rk4 reiterations
    for step in range(int(L / dx)):
        [new_T, new_Tprime] = rk4(T[-1], Tprime[-1])
        T.append(new_T)
        Tprime.append(new_Tprime)

    # collect value of T(L)
    T_LList.append(T[-1] - T_L)

plt.figure()
plt.plot(sList, T_LList)
# plt.legend(['0.2','0.05'])
plt.ylabel('Boundary condition error [K]')
plt.xlabel('Value s set in T\'(0) = s')
plt.title('Approximation of T\'(0)')

minimumError = 1000
index = 0
indexCount = 0
for num in T_LList:
    if (abs(num) < abs(minimumError)):
        minimumError = num
        index = indexCount
    indexCount = indexCount + 1

print("Closest root: " + str(minimumError) + " at s = " + str(sList[index]))
s2 = sList[index]

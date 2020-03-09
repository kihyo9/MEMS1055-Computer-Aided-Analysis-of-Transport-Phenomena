# Shooting Method + rk4 + plots
import matplotlib.pyplot as plt
import numpy as np

# constants
k = 1
L = 10
sigma = 2.7e-9
h = 0.05
T_inf = 200
T_0 = 300
T_L = 400
C = (-1 / k) * (h * T_inf + sigma * T_inf ** 4)  # -14.32

# numerical solution parameters
dx = 0.1
leftEnd, rightEnd = -250, -39
sList = np.linspace(leftEnd, rightEnd, (rightEnd - leftEnd)*10 + 1)  # Guesses for dt/dx @ x=0
T_LList = []  # List of errors from the exact BC: T(L) = 400


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

#shooting method - calculate T(L) values for various s
for s in sList:
    s = round(s,2)
    print(s)
    # initial conditions
    T = [T_0]
    Tprime = [s]  # dT/dx

    # rk4 iterations
    for step in range(int(L / dx)):
        [new_T, new_Tprime] = rk4(T[-1], Tprime[-1])
        T.append(new_T)
        Tprime.append(new_Tprime)

    # collect value of T(L)
    T_LList.append(T[-1] - T_L)

#find and print roots
roots = []
for i, T in enumerate(T_LList):
    if i > 0:
        if T_LList[i] * T_LList[i-1] < 0:
            if abs(T_LList[i]) < abs(T_LList[i-1]):
                roots.append([i, round(sList[i],2)])
            else:
                roots.append([i-1,round(sList[i-1],2)])

for i, root in enumerate(roots):
    print('Root ' + str(i) + ': ' + str(root[1]))

#plot s vs error for T(L)
plt.figure()
markers_on = [root[0] for root in roots]
plt.plot(sList, T_LList,marker='o',markevery=markers_on)
# plt.legend(['0.2','0.05'])
plt.ylabel('Boundary condition T(L) error [K]')
plt.xlabel('s')
plt.title('Approximation of T\'(0) = s')

#plotting solutions for roots of s plot
for s in [root[1] for root in roots]:
    print(s)
    # initial conditions
    T = [T_0]
    Tprime = [s]  # dT/dx

    # rk4 reiterations
    for step in range(int(L / dx)):
        [new_T, new_Tprime] = rk4(T[-1], Tprime[-1])
        T.append(new_T)
        Tprime.append(new_Tprime)

    plt.figure()
    plt.plot(np.linspace(0,L,int((L/dx)+1)), T)
    # plt.legend(['0.2','0.05'])
    plt.ylabel('Temperature [K]')
    plt.xlabel('Position [m]')
    plt.title('Approximation of T(x) with s = ' + str(s))

plt.show()
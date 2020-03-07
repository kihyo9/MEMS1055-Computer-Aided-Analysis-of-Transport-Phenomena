import matplotlib.pyplot as plt
import numpy as np

#rk4 - second root as s

k = 1
L = 10
sigma = 2.7e-9
h=0.05
T_inf = 200
T_0 = 300
T_L = 400
C = (-1/k)*(h*T_inf + sigma*T_inf**4)

dx = 0.0001
s = s2 # = dt/dx @ x=0

#define rk4 calculation
def rk4(current_T, current_Tprime):
  T1 = current_Tprime*dx
  Tprime1 = Tprimeprime(current_T)*dx

  T2 = (current_Tprime + 0.5*Tprime1)*dx
  Tprime2 = Tprimeprime(current_T + 0.5*T1)*dx

  T3 = (current_Tprime + 0.5*Tprime2)*dx
  Tprime3 = Tprimeprime(current_T + 0.5*T2)*dx

  T4 = (current_Tprime + Tprime3)*dx
  Tprime4 = Tprimeprime(current_T + T3)*dx

  new_T = current_T + (1/6)*(T1+2*T2+2*T3+T4)
  new_Tprime = current_Tprime + (1/6)*(Tprime1+2*Tprime2+2*Tprime3+Tprime4)

  return new_T, new_Tprime

#define d2T/dx2
def Tprimeprime(current_T):
  try:
    return (h/k)*current_T + (sigma/k)*current_T**4 + C
  except:
    print(str(current_T))
    raise

#initial conditions
T = [T_0]
Tprime = [s] #dT/dx

#reiterating rk4
for step in range(int(L/dx)):
  [new_T, new_Tprime] = rk4(T[-1], Tprime[-1])
  T.append(new_T)
  Tprime.append(new_Tprime)

plt.figure()
plt.plot(np.linspace(0,L,(L/dx)+1), T)
# plt.legend(['0.2','0.05'])
plt.ylabel('Temperature [K]')
plt.xlabel('Position [m]')
plt.title('Approximation of T(x) with s = ' + str(s2))

print(str(T[-1]))

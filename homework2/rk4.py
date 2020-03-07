import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

###############################################################
# RK4
# constants
l = 0.6
g = -9.81
theta_start = 0.1

def rk4(theta, theta_dot):
    current_theta = theta[-1]
    current_omega = theta_dot[-1]

    theta1 = current_omega * dt
    omega1 = dOmega(current_theta) * dt

    theta2 = (current_omega + omega1*0.5) * dt
    omega2 = dOmega(current_theta + theta1*0.5) * dt

    theta3 = (current_omega + omega2*0.5) * dt
    omega3 = dOmega(current_theta + theta2*0.5) * dt

    theta4 = (current_omega + omega3) * dt
    omega4 = dOmega(current_theta + theta3) * dt

    theta.append(current_theta + (1/6)*(theta1 + 2*theta2 + 2*theta3 + theta4))
    theta_dot.append(current_omega + (1/6)*(omega1 + 2*omega2 + 2*omega3 + omega4))

def dOmega(theta_dd):
    return (g/l)*np.sin(theta_dd)

# dts = np.linspace(0.7, 0.01, 70)
dts = np.linspace(0.04, 0.001, 40)
sums = []

fig1 = plt.figure()
for dt in dts:
    dt = round(dt,3)
    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]

    for x in range(int(10. / dt)):
        x_dt.append(dt * x)
        rk4(theta, theta_dot)

    print(dt)
    sum = 0
    analytic = hp.analyticalSol(x_dt,theta_start)
    for i in range(len(theta)):
        sum = sum + (theta[i] - analytic[i]) ** 2
    sums.append((1 / len(theta)) * sum)

    # plotting
    # if dt in [0.4, 0.25]:
    if dt in [0.007]:
        plt.plot(x_dt, theta)

plt.title('Solutions of various dt')
plt.plot(np.linspace(0,10,1001), hp.analyticalSol(np.linspace(0,10,1001),theta_start))
plt.legend(['0.007[s]','Analytic'],loc='lower left')
# plt.legend(['0.4[s]','0.25[s]','Analytic'],loc='lower left')
plt.ylabel('Amplititude')
plt.xlabel('time [s]')


plt.figure()
plt.title('Timestep vs Residuals')
plt.plot(dts, sums, marker='o')
plt.xlabel('Timestep [s]')
plt.ylabel('Residual value')
plt.show()
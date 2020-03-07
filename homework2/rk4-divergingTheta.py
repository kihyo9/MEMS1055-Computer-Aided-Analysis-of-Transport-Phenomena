import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

# constants
l = 0.6
g = -9.81
dt = 0.02
theta_starts = np.linspace(0.01, np.pi, 314)

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
    if isSin:
       return (g/l)*np.sin(theta_dd)
    else:
        return (g / l) * theta_dd

#####################################################################
# rk4


difference = []
isSin = True

for theta_start in theta_starts:
    dt = round(dt,3)
    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]
    isSin = True

    for x in range(int(10. / dt)):
        x_dt.append(dt * x)
        rk4(theta, theta_dot)

    equation1theta = theta.copy()

    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]
    isSin = False
    for x in range(int(10. / dt)):
        x_dt.append(dt * x)
        rk4(theta, theta_dot)

    equation4theta = theta.copy()

    print(theta_start)

    diff = 0
    for eq1, eq4 in zip(equation1theta, equation4theta):
        diff += (eq1 - eq4)**2

    difference.append(diff/len(equation1theta))

plt.plot(theta_starts, difference)
plt.title('Starting angle vs Solution differences')
# plt.legend(['0.01[rad]','0.1[rad]','0.5[rad]','1[rad]'])
plt.ylabel('Solution differences')
plt.xlabel('Starting angle')
plt.show()
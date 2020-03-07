import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

# constants
l = 0.6
g = -9.81
dt = 0.02
theta_starts = np.linspace(0.01, 0.5, 100)

#####################################################################
# explicit


difference = []

for theta_start in theta_starts:
    dt = round(dt,3)
    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]
    for x in range(int(10. / dt)):
        x_dt.append(dt * x)

        next_theta = theta[-1] + theta_dot[-1] * dt
        next_theta_dot = theta_dot[-1] + (g / l) * np.sin(theta[-1]) * dt
        theta.append(next_theta)
        theta_dot.append(next_theta_dot)

    equation1theta = theta.copy()

    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]
    for x in range(int(10. / dt)):
        x_dt.append(dt * x)

        next_theta = theta[-1] + theta_dot[-1] * dt
        next_theta_dot = theta_dot[-1] + (g / l) * theta[-1] * dt
        theta.append(next_theta)
        theta_dot.append(next_theta_dot)

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
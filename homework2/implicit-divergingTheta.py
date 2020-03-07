import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

# constants
l = 0.6
g = -9.81
dt = 0.02
theta_starts = np.linspace(0.01, 1, 100)

#####################################################################
# implicit


difference = []

for theta_start in theta_starts:
    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]

    for x in range(int(10. / dt)):
        x_dt.append(dt * x)

        # The system of equations to solve
        # next_theta = theta + next_theta_dot * dt
        # next_theta_dot = theta_dot + (g / l) * np.sin(next_theta) * dt

        theta_func = lambda x: theta[-1] + x[1] * dt
        theta_dot_func = lambda x: theta_dot[-1] + (g/l) * x[0] * dt
        soe = [theta_func, theta_dot_func]

        dtheta_func = lambda x: 2*(x[0] - theta_func(x)) + 2*(x[1] - theta_dot_func(x))*(-1*(g/l) * dt)
        dtheta_func_dot = lambda x: 2*(x[0] - theta_func(x))*(-1*dt) + 2*(x[1] - theta_dot_func(x))
        partials = [dtheta_func, dtheta_func_dot]
        [next_theta, next_theta_dot] = hp.gradientDescent(soe, [theta[-1], theta_dot[-1]], hp.errorDef1, partials, 0.01, 0.0000001)

        theta_dot.append(next_theta_dot)
        theta.append(next_theta)

    equation1theta = theta.copy()

    theta = [theta_start]
    theta_dot = [0]
    x_dt = [0]

    for x in range(int(10. / dt)):
        x_dt.append(dt * x)

        # The system of equations to solve
        # next_theta = theta + next_theta_dot * dt
        # next_theta_dot = theta_dot + (g / l) * np.sin(next_theta) * dt

        theta_func = lambda x: theta[-1] + x[1] * dt
        theta_dot_func = lambda x: theta_dot[-1] + (g/l) * np.sin(x[0]) * dt
        soe = [theta_func, theta_dot_func]

        dtheta_func = lambda x: 2*(x[0] - theta_func(x)) + 2*(x[1] - theta_dot_func(x))*(-1*(g/l) * np.cos(x[0]) * dt)
        dtheta_func_dot = lambda x: 2*(x[0] - theta_func(x))*(-1*dt) + 2*(x[1] - theta_dot_func(x))
        partials = [dtheta_func, dtheta_func_dot]
        [next_theta, next_theta_dot] = hp.gradientDescent(soe, [theta[-1], theta_dot[-1]], hp.errorDef1, partials, 0.01, 0.0000001)

        theta_dot.append(next_theta_dot)
        theta.append(next_theta)

    equation4theta = theta.copy()

    print(theta_start)

    diff = 0
    for eq1, eq4 in zip(equation1theta, equation4theta):
        diff += (eq1 - eq4)**2

    difference.append(diff/len(equation1theta))

plt.plot(theta_starts, difference,marker='o')
plt.title('Starting angle vs Solution differences')
# plt.legend(['0.01[rad]','0.1[rad]','0.5[rad]','1[rad]'])
plt.ylabel('Solution differences')
plt.xlabel('Starting angle')
plt.show()

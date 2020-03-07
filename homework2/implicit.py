import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

# constants
l = 0.6
g = -9.81
theta_start = 0.1
#####################################################################
# implicit residuals
dts = np.linspace(0.50, 0.01, 50)
# dts = np.linspace(0.02, 0.001, 20)

sums = []

for dt in dts:
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

    print(dt)
    sum = 0
    analytic = hp.analyticalSol(x_dt, theta_start)
    for i in range(len(theta)):
        ratio = dt / 0.01
        sum = sum + (theta[i] - analytic[i]) ** 2
    sums.append((1 / len(theta)) * sum)

    # plotting
    # if dt == 0.001 or dt == 0.008:
    if dt == 0.02 or dt == 0.005:
        plt.plot(x_dt, theta)

plt.title('Solutions of various dt')
plt.plot(np.linspace(0,10,1001), hp.analyticalSol(np.linspace(0,10,1001), theta_start))
plt.legend(['0.02[s]','0.005[s]','Analytic'])
plt.ylabel('Amplititude')
plt.xlabel('time [s]')

plt.figure()
plt.title('Timestep vs Residuals')
plt.plot(dts, sums, marker='o')
plt.xlabel('Timestep [s]')
plt.ylabel('Residual value')
plt.show()
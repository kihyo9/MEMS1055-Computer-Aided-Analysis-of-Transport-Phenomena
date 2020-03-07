import matplotlib.pyplot as plt
import numpy as np
import helpers as hp

# constants
l = 0.6
g = -9.81
theta_start = 0.1

#####################################################################
# explicit residuals
dts = np.linspace(0.49, 0.01, 49)
# dts = np.linspace(0.04, 0.001, 40)
res = []
sums = []

fig1 = plt.figure()
for dt in dts:
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

    print(dt)
    sum = 0
    analytic = hp.analyticalSol(x_dt,theta_start)
    for i in range(len(theta)):
        sum = sum + (theta[i] - analytic[i]) ** 2
    sums.append((1 / len(theta)) * sum)

    # plotting
    if dt == 0.01 or dt == 0.03:
    # if dt in [0.007,0.005, 0.002]:
        plt.plot(x_dt, theta)

plt.title('Solutions of various dt')
plt.plot(np.linspace(0,10,1001), hp.analyticalSol(np.linspace(0,10,1001),theta_start))
# plt.legend(['0.007[s]','0.005[s]','0.002[s]','Analytic'])
plt.legend(['0.03[s]','0.01[s]','Analytic'])
plt.ylabel('Amplititude')
plt.xlabel('Time [s]')

plt.figure()
plt.title('Timestep vs Residuals')
plt.plot(dts, sums, marker='o')
plt.xlabel('Timestep [s]')
plt.ylabel('Residual value')
plt.show()
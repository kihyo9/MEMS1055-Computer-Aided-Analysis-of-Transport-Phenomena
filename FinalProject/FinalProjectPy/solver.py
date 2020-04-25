import matplotlib.pyplot as plt
import time
from datetime import datetime


def plotting(fileName, dx_data, dt_data, gamma):

    f = open(fileName, "r")
    fileContent = f.readlines()
    f.close()

    plt.figure()
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    legend = []
    plt.title("Gamma = %.2f" % gamma)

    l = len(fileContent)
    for i in range(len(fileContent)):
        timestepPlot(i, fileContent, dx_data, dt_data, legend, l)
    # timestepPlotAll(fileContent, dx_data, dt_data, legend, l)

    plt.legend(legend,loc='center left',bbox_to_anchor=(1, 0.5))
    plt.xlabel("x [m]")
    plt.ylabel("c")
    plt.show()

def plotting2(fileName, title):

    f = open(fileName, "r")
    fileContent = f.readlines()
    f.close()

    xMin, xMax, dx = map(float, fileContent[0].split())
    tMin, tMax, dt = map(float, fileContent[1].split())
    dx_steps = len(fileContent[3].split())
    dt_steps = int(len(fileContent)/2) - 1
    dx_data = [xMin + dx * i for i in range(dx_steps)]
    dt_data = [tMin + dt * i for i in range(dt_steps)]

    plt.figure()
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    legend = []
    plt.title(title)

    l = len(fileContent)
    for i in range(dt_steps):
        timestepPlot(2 + 2*i, fileContent, dx_data, dt_data, legend, l)

    plt.legend(legend,loc='center left',bbox_to_anchor=(1, 0.5))
    plt.xlabel("x [m]")
    plt.ylabel("c")
    plt.show()

def printStability(dx, dt, u, gamma):
    print("\n_____Starting Simulation_____\n")
    print("Gamma = {}, u = {}, timestep size = {}".format(gamma, u, dt))
    r = gamma*dt/dx**2
    print("CFL is %f (must be < %f)" % (dt*u/dx, (1 - 4*r**2)**0.5))
    print("r is %f (must be less than 0.5)" % r)
    print("G is %f (must be greater than 0)" % (2*gamma/dt - u**2))

def timestepPlot(stepline, fileContent, dx_data, dt_data, legend, l):
    timestep, time = fileContent[stepline].split()
    fdata = fileContent[stepline+1].split()
    floatdata = [float(s) for s in fdata]
    plt.plot(dx_data, floatdata)
    index = sequenceGenerator(stepline, len(dt_data)-1)
    legend.append("%s: %s" %(timestep,time))

def timestepPlotAll(fileContent, dx_data, dt_data, legend, l):
    startplot = 0
    timestepPlot(startplot, fileContent, dx_data, dt_data, legend, l)
    startplot = 1
    while startplot < len(dt_data):
        timestepPlot(startplot, fileContent, dx_data, dt_data, legend, l)
        startplot = int(round(startplot*2.16,0))
    else:
        startplot = -1
        timestepPlot(startplot, fileContent, dx_data, dt_data, legend, l)

def sequenceGenerator(i, max):
    # if i number of lines are written to file, then write a line again when n = i-th term sequence
    if i == 0:
        return 0
    elif i == 1:
        return 1
    else:
        ans = 1
        for _ in range(i-1):
            ans = int(round(ans*2.16,0))
        if ans > max:
            return max
        else:
            return ans

def writeToFile(fileName, t_step, n, linesWritten, dt_steps):
    if sequenceGenerator(linesWritten, dt_steps) == n:
        f = open(fileName, "a")
        for i,num in enumerate(t_step):
            if i == len(t_step) - 1:
                f.write(str(num) + "\n")
            else:
                f.write(str(num) + ", ")
        f.close()
        return True
    else:
        return False

def performanceTime(starttime, solverType, gamma, dt_steps,u,L):
    g = open("performance.txt","a")
    g.write("Start time: {}\n".format(datetime.now()))
    g.write("{}: gamma = {}, timesteps = {}, u = {}, L = {}\n".format(solverType,gamma,dt_steps,u,L))
    g.write(str('{}s\n\n'.format(time.time() - starttime)))
    g.close()

if __name__ == "__main__":
    f = open("sequence.txt","w")
    for i in range(50):
        f.write("%d\n" % (sequenceGenerator(i, 1e7)))

class cat:
    def yeeet(self, yeet):
        print(yeet)
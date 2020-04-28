import matplotlib.pyplot as plt
import numpy as np
import implicit
import explicit

k = 5. # doesnt matter
R = 0.25
dr = 0.005
dr_steps = int(round(R/dr, 0))
BiList = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 20, 50]
# BiList = [1, 2, 5, 20, 50]
T_inf = 300.
T_0 = 400.
rho = 7300.
cp = 460.
alpha = 23e-6  # of cylinder #k_fluid/(cp*rho)
G = 0.20  # set this
dt = G * dr ** 2 / alpha

# explicit
figure511points_ex = []
plt.figure()
plt.title("Explicit")
Q_0 = rho*cp*R**2 *np.pi*(T_0 - T_inf) # per length of cylinder
Bi2Fo = lambda h,alpha,t,k: (h/k)**2 * alpha*t
for Bi in BiList:
    print("Bi: " + str(Bi))
    h = k/(Bi*R)

    explicit.explicitSolver(k,h,"figure512-explicit-data.txt")
    f = open("figure512-explicit-data.txt", "r")
    fileContentf = f.readlines()

    Qs = []
    Bi2Fos = []

    for i,line in enumerate(fileContentf):
        if i % 50 == 0:
            lineData = line.split(", ")
            Q = 0
            for j,gridpoint in enumerate(lineData):
                if j == 0:
                    continue
                Q += rho*cp*np.pi*(T_0 - float(gridpoint))*((2*j-1)*dr**2)
            Qs.append(Q/Q_0)

            Bi2Fos.append(Bi2Fo(h,alpha,dt*i,k))

    g = open("debug.txt","a")
    g.write(str(Bi) + "\n")
    g.write(str(Qs) + "\n")
    g.write(str(Bi2Fos) + "\n\n")
    g.close()

    f.close()
    plt.plot(Bi2Fos, Qs)

plt.xlabel("Bi^2*Fo")
plt.ylabel("Q/Q_0")
plt.xscale("log")
plt.legend(["Bi = " + str(Bi) for Bi in BiList])

# # implicit
# figure511points_im = []
# plt.figure()
# plt.title("Implicit")
# for iBi in inverseBi:
#     h = k/(iBi*R)
#     print("1/Bi: " + str(k / (h * R)))
#     implicit.implicitSolver(k, h,"figure511-implicit-data.txt")
#     g = open("figure511-implicit-data.txt", "r")
#     fileContentg = g.readlines()
#     timestep_chosen = int(len(fileContentg)/2)
#     gdata = fileContentg[timestep_chosen].split(", ")
#     figure511points_im.append([(float(gdata[int(round(radius*dr_steps,0))]) - T_inf)/(float(gdata[0]) - T_inf) for radius in radii])
#     g.close()
# figure511points_imT = list(zip(*figure511points_im))
# for i,radius in enumerate(radii):
#     plt.plot(inverseBi,figure511points_imT[i])
# plt.xlabel("1/Bi")
# plt.ylabel("theta/theta_0")
# plt.xscale("log")
# plt.legend(["r/R = " + str(radius) for radius in radii])

plt.show()
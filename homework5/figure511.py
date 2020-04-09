import matplotlib.pyplot as plt
import numpy as np
import implicit
import explicit

k = 5. # of cylinder
R = 0.25
dr = 0.005
dr_steps = int(round(R/dr, 0))
radii = [0.2,0.4,0.6,0.8,0.9,1.]
inverseBi = [0.02, 0.1, 0.5, 1, 5, 20, 50, 100, 200, 500, 1000]
T_inf = 300

# explicit
figure511points_ex = []
plt.figure()
plt.title("Explicit")
for iBi in inverseBi:
    h = k/(iBi*R)
    print("1/Bi: " + str(k / (h * R)))
    explicit.explicitSolver(k,h,"figure511-explicit-data.txt")
    f = open("figure511-explicit-data.txt", "r")
    fileContentf = f.readlines()
    timestep_chosen = int(len(fileContentf)/2)
    fdata = fileContentf[timestep_chosen].split(", ")
    figure511points_ex.append([(float(fdata[int(round(radius*dr_steps,0))]) - T_inf)/(float(fdata[0]) - T_inf) for radius in radii])
    f.close()
figure511points_exT = list(zip(*figure511points_ex))
for i,radius in enumerate(radii):
    plt.plot(inverseBi,figure511points_exT[i])
plt.xlabel("1/Bi")
plt.ylabel("theta/theta_0")
plt.xscale("log")
plt.legend(["r/R = " + str(radius) for radius in radii])

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
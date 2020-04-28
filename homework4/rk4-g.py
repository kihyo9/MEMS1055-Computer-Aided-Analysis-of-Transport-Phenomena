import matplotlib.pyplot as plt
import numpy as np
print("right")
### Functions for rk4 ###
def rk4(f,G,H):
    current_F = f
    current_G = G[-1]
    current_H = H[-1]

    G1 = current_H * dt
    H1 = dH(current_F,current_H) * dt

    G2 = (current_H + H1*0.5) * dt
    H2 = dH(current_F, current_H + H1*0.5) * dt

    G3 = (current_H + H2*0.5) * dt
    H3 = dH(current_F, current_H + H2*0.5) * dt

    G4 = (current_H + H3) * dt
    H4 = dH(current_F, current_H + H3) * dt

    G.append(current_G + (1/6)*(G1 + 2*G2 + 2*G3 + G4))
    H.append(current_H + (1/6)*(H1 + 2*H2 + 2*H3 + H4))

def dH(f,h):
    return -(Pr)*f*h


### rk4 + shooting method ###
# Find the root(s) of the plot of different values of s
# The diffeq is : h' + Pr*f*h=0
#               : g' = h
# BC1: g(0) = 0
# BC2: g(inf) = 0

dt = 0.01
sList = np.linspace(0,1,201)
F = []
PrList = [0.5, 1, 2]

# Get the F values from file
f = open("Fplot.txt", "r")
fileContent = f.readlines()
for line in fileContent:
    F.append(float(line))

roots = []
for prandtl in PrList:
    Pr = prandtl
    G_inf = []
    for s in sList:
        G = [0]
        H = [s]
        t_axis = [0]
        for i,x in enumerate(range(int(10. / dt))):
            dt = round(dt, 3)
            t_axis.append(dt * x)
            rk4(F[i],G,H)
        G_inf.append(G[-1]-1)

    ### find the roots of the above plot ###
    first = True
    last = 0
    for num, el in enumerate(G_inf):
        if first:
            last = el
            first = False
            continue
        if el*last < 0 or el == 0:
            roots.append([sList[num],num])
        last = el

### for what eta does F'(eta) = 0.99*F'(inf)? ###
def magicEta(fG,ft_axis):
    first = True
    last = 0
    for num,el in enumerate(fG):
        if first:
            last = el
            first = False
            continue
        if el > 0.99:
            return ft_axis[num]
        last = el
    return None

### plot the solution with the estimate for s ###
G = []
t_axis = []
for num,s in enumerate(roots):
    plt.figure()
    legend = []
    Pr = PrList[num]
    t_axis = [0]
    G = [0]
    H = [s[0]]
    for i,x in enumerate(range(int(10. / dt))):
        dt = round(dt, 3)
        t_axis.append(dt * x)
        try:
            rk4(F[i],G,H)
        except:
            print("error: " + str(i))
    plt.plot(t_axis,G)
    plt.plot(t_axis,H)
    legend.append('G')
    legend.append('dG/dEta')
    plt.title('Plots of G and dG/dEta; Pr = ' + str(Pr))
    plt.legend(legend)
    plt.ylabel('value')
    plt.xlabel('eta')
    print("Pr = " + str(Pr) + ", magic eta  = " + str(magicEta(G, t_axis)))
    print("s = " + str(roots[num][0]) + "\n")

print("s = " + str(roots[0]))

# show plot
plt.show()
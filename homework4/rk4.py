import matplotlib.pyplot as plt
import numpy as np


### Functions for rk4 ###
def rk4(F,G,H):
    current_F = F[-1]
    current_G = G[-1]
    current_H = H[-1]

    F1 = current_G * dt
    G1 = current_H * dt
    H1 = dH(current_F,current_H) * dt

    F2 = (current_G + G1*0.5) * dt
    G2 = (current_H + H1*0.5) * dt
    H2 = dH(current_F + F1*0.5, current_H + H1*0.5) * dt

    F3 = (current_G + G2*0.5) * dt
    G3 = (current_H + H2*0.5) * dt
    H3 = dH(current_F + F2*0.5, current_H + H2*0.5) * dt

    F4 = (current_G + G3) * dt
    G4 = (current_H + H3) * dt
    H4 = dH(current_F + F3, current_H + H3) * dt

    F.append(current_F + (1/6)*(F1 + 2*F2 + 2*F3 + F4))
    G.append(current_G + (1/6)*(G1 + 2*G2 + 2*G3 + G4))
    H.append(current_H + (1/6)*(H1 + 2*H2 + 2*H3 + H4))

def dH(f,h):
    return -f*h/2

### rk4 + shooting method ###
# Find the root(s) of the plot of different values of s

difference = []
dt = 0.01
t_axis = [0]
sList = np.linspace(0,5,201)
G_inf = []

for s in sList:
    F = [0]
    G = [0]
    H = [s]
    for x in range(int(10. / dt)):
        dt = round(dt, 3)
        t_axis.append(dt * x)
        rk4(F,G,H)
    G_inf.append(G[-1] - 1)

plt.figure()
plt.plot(sList, G_inf)
plt.title("Plot of guesses for f'(inf) = s")
plt.ylabel('G_inf - 1')
plt.xlabel('guess for s')

### find the roots of the above plot ###
roots = []
first = True
last = 0
for num, el in enumerate(G_inf):
    if first:
        last = el
        first = False
        continue
    if el*last < 0:
        roots.append([sList[num],num])
    last = el

### plot the solution with the estimate for s ###
plt.figure()
legend = []
F = []
G = []
t_axis = []
for num,s in enumerate(roots):
    t_axis = [0]
    F = [0]
    G = [0]
    H = [s[0]]
    for x in range(int(10. / dt)):
        dt = round(dt, 3)
        t_axis.append(dt * x)
        rk4(F,G,H)
    plt.plot(t_axis,F)
    plt.plot(t_axis,G)
    plt.plot(t_axis,H)
    legend.append('F')
    legend.append('dF/dEta')
    legend.append('d2F/dEta2')

plt.title('Plots of F,G and H')
plt.legend(legend)
plt.ylabel('value')
plt.xlabel('eta')

# shear stress

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

print("s = " + str(roots[0]))
print("magic eta  = " + str(magicEta(G,t_axis)))

### save F value to file ###
f = open("Fplot.txt", "w")
for i in F:
    f.write(str(i) + "\n")
f.close()

### show plots ###
plt.show()
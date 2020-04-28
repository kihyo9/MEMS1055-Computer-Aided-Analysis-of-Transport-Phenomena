import numpy as np
import matplotlib.pyplot as plt

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver

def thomasAlgorithmSolver(matrix):
    a = matrix[0]
    b = matrix[1]
    c = matrix[2]
    d = matrix[3]

    n = len(d)

    cc = []
    dd = []
    for i in range(n):
        if i == 0:
           cc.append(c[i]/b[i])
           dd.append(d[i] / b[i])
        elif i == n - 1:
            dd.append((d[i] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1]))
        else:
            cc.append(c[i]/(b[i] - a[i]*cc[i-1]))
            dd.append((d[i] - a[i] * dd[i - 1]) / (b[i] - a[i] * cc[i - 1]))

    x = [0 for i in range(n)]
    for i in range(n):
        q = n-1-i

        if q == n-1:
            x[q] = dd[q]
        else:
            x[q] = dd[q] - cc[q]*x[q+1]
    return x

if __name__ == "__main__":
    f = open("Tdata-explicit.txt", "r")
    g = open("Tdata-implicit.txt", "r")

    fileContentf = f.readlines()
    fileContentg = g.readlines()

    differences = []
    for linef,lineg in zip(fileContentf,fileContentg):
        fdata = linef.split(", ")
        gdata = lineg.split(", ")

        dif = 0
        for fpoint, gpoint in zip(fdata, gdata):
            dif += abs(float(fpoint) - float(gpoint))
        differences.append(dif/len(fdata))

    f.close()
    g.close()

    plt.figure()
    plt.plot([i*1e-5 for i in range(len(differences))], differences)
    plt.title('Average difference between implicit and explicit solutions per gridpoint')
    plt.ylabel('Temperature [K]')
    plt.xlabel('Time [s]')
    plt.show()


    '''
    A = [0.,1.,1.,1.,1.]
    B = [2.,2.,2.,2.,2.]
    C = [3.,3.,3.,3.,0.]
    D = [7.,7.,7.,7.,7.]
    matrix = [A,B,C,D]
    print(thomasAlgorithmSolver(matrix))
    aa = np.array([[2.,3.,0.,0., 0,],[1.,2.,3.,0., 0.],[0.,1.,2.,3.,0],[0.,0.,1.,2.,3.],[0.,0.,0.,1.,2.]])
    bb = np.array([7.,7.,7.,7.,7.])
    print(np.linalg.solve(aa,bb))
    '''
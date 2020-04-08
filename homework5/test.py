import numpy as np


## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    nf = len(d)  # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))  # copy arrays
    for it in range(1, nf):
        mc = ac[it - 1] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        dc[it] = dc[it] - mc * dc[it - 1]

    xc = bc
    xc[-1] = dc[-1] / bc[-1]

    for il in range(nf - 2, -1, -1):
        xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

    return xc

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
def thomasAlgorithmSolver2(matrix):
    a = matrix[0]
    b = matrix[1]
    c = matrix[2]
    d = matrix[3]

    n = len(d)

    w = [0 for i in range(n)]
    x = [0 for i in range(n)]
    for i in range(n):
        if i == 0:
            continue
        else:
            w.append(a[i]/b[i-1])
            b[i] -= w[i]*c[i-1]
            d[i] -= w[i]*d[i-1]
            a[i] = 0

    for i in range(n):
        if i == 0:
            x[n-1] = d[n-1]/b[n-1]
        else:
            x[n-1-i] = (d[n-1-i] - c[n-1-i]*x[n-i])/b[n-1-i]

    return x
A = [0.,1.,1.,1.]
B = [2.,2.,2.,2.]
C = [3.,3.,3.,0.]
D = [7.,7.,7.,7.]
matrix = [A,B,C,D]
print(thomasAlgorithmSolver(matrix))
print(TDMAsolver(matrix[0],matrix[1],matrix[2],matrix[3]))
print(thomasAlgorithmSolver2(matrix))
aa = np.array([[2.,3.,0.,0.],[1.,2.,3.,0.],[0.,1.,2.,3.],[0.,0.,1.,2.]])
bb = np.array([7.,7.,7.,7.])
print(np.linalg.solve(aa,bb))
# answer: 2.5455    1.9091   -4.4545   10.1818
# current output: [0.9166666666666666, 0.4583333333333333, 0.0, 2.5454545454545454]
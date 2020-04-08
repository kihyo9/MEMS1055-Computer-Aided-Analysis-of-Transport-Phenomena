def thomasAlgorithmSolver(matrix):
    a = matrix[0]
    b = matrix[1]
    c = matrix[2]
    d = matrix[3]

    n = len(d)
    cc = []
    dd = []
    x = []

    # cc
    for i in range(n-1):
        if i == 0:
            cc.append(c[i]/b[i])
        else:
            cc.append(c[i] / (b[i] - a[i]*cc[i-1]))

    # dd
    for i in range(n):
        if i == 0:
            dd.append(d[i]/b[i])
        else:
            dd.append((d[i]-a[i]*dd[i-1])/(b[i]-a[i]*cc[i-1]))

    # x
    for i in range(n):
        i = (n - 1) - i

        if i == n - 1:
            x.append(dd[i])
        else:
            x = [dd[i]/(cc[i]*x[-1])] + x

    return x
matrix = [[1.,1.,1.,1.],[2.,2.,2.,2.],[3.,3.,3.,3.],[7.,7.,7.,7.]]
print(thomasAlgorithmSolver(matrix))
# answer: 2.5455    1.9091   -4.4545   10.1818
# current output: [0.9166666666666666, 0.4583333333333333, 0.0, 2.5454545454545454]
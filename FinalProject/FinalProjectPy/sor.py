import numpy as np
import time

def sor(matrix):
    pass



size = 10
start = time.time()
rana = np.random.rand(size,size)
ranb = np.random.rand(size)
x = np.linalg.solve(rana,ranb)
print(x)
finish = time.time()
print(finish - start)


# a = np.array([[4,-1,-6,0],[-5,-4,10,8],[0,9,4,-2],[1,0,-7,5]])
# b = np.array([2,21,-12,-6])
# x = np.linalg.solve(a,b)
# print(x)

start = time.time()
w = 1.5
answer = [0 for _ in range(size)]
count = 0
while(abs(np.sum(rana.dot(answer) - ranb)) > 0.0001):
    count +=1
    for i,array in enumerate(rana):
        sig = 0
        for j in range(size):
            if i == j:
                continue
            else:
                sig += array[j]*answer[j]
        answer[i] = (1-w)*answer[i] + (w/array[i]) * (ranb[i] - sig)
    print(str(count), ", ", answer)
finish = time.time()
print(answer)
print(finish - start)

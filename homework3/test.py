# temperature dist of the bar
# shooting+rk4 is finished
# finite difference with thomas method
# determine optimum dx for each case
## maybe find the difference in solution with each successive decrease in dx?

import numpy as np
import matplotlib.pyplot as plt

def error(x,s):
    y = []
    for q in x:
        y.append((2*q+3*s - 4)**2 + (q+2*s+3*s - 4)**2)
    return y

plt.figure()
arr = np.linspace(-11,9,201)
plt.plot(arr,error(list(arr),3))
plt.plot(arr,error(list(arr),2))
plt.plot(arr,error(list(arr),1))
plt.plot(arr,error(list(arr),0))
plt.plot(arr,error(list(arr),-1))
plt.plot(arr,error(list(arr),-2))
plt.legend(['3','2','1','0','-1','-2'])
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from sys import exit


rho = np.loadtxt("rhoz.dat")
x = rho[:,0]
y = rho[:,1] 

dx = x[1] - x[0]
print( dx * sum(y))

plt.plot(x,y)
plt.show()

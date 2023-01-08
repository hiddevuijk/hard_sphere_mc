import numpy as np
import matplotlib.pyplot as plt
from sys import exit


rho = np.loadtxt("rhoz.dat")
x = rho[:,0]
y = rho[:,1] / (1200/2.0)

plt.plot(x,y)
plt.show()

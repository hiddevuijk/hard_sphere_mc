import numpy as np
import matplotlib.pyplot as plt
from sys import exit


gr = np.loadtxt("gr.dat")
r = gr[:,0]
y = gr[:,1]

plt.plot(r,y, color="blue", label="mc")


gr = np.loadtxt("../hard_sphere_ed_bd/data/gr.dat")
r = gr[:,0]
y = gr[:,1] 

plt.plot(r,y, color="red", label="ed_bd")

plt.legend()
plt.show()

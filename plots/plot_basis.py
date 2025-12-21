import numpy as np
import matplotlib.pyplot as plt
import csv

data = np.loadtxt("data/basis.csv", delimiter=",")
u = data[:,0]
vals = data[:,1:]

for i in range(vals.shape[1]):
    plt.plot(u, vals[:,i], label=f"N{i}")

plt.legend()
plt.show()

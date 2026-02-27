import numpy as np
import matplotlib.pyplot as plt

pressData = np.loadtxt('ld_2D_NonInter_pressure.dat')
t = pressData[:, 0]
p = pressData[:, 1]


N = 1000
kT = 1.0
L = 10.0
p_ideal = N * kT / (L*L)

plt.axhline(p_ideal, linestyle='--', color='black')

plt.plot(t, p)
plt.show()
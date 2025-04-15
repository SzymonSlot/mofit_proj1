import numpy as np
import matplotlib.pyplot as plt
from funkcje import *


# zadanie 1 Verlet

dttab = [ 0.02, 0.016, 0.01, 0.008, 0.005, 0.002, 1e-3]
m=1
x10 = np.zeros(len(dttab))
t_list = []
x_list = []
v_list = []
E_list = []
i=0


for dt in dttab:
    t = np.arange(0, 100, dt)
    x = np.zeros(len(t))
    v = np.zeros(len(t))
    a = np.zeros(len(t))

    x[0] = 2.8
    v[0] = 0

    for n in range(0, len(t)-1):
        x[n+1], v[n+1] = verlet(x[n], v[n], dt, m)


    E = m*v**2/2 + V(x)

    t_list.append(t)
    x_list.append(x)
    v_list.append(v)
    E_list.append(E)

    idx = np.abs(t - 10).argmin()
    x10[i] = x[idx]
    i+=1

# plotting(t_list,x_list,v_list,E_list,dttab,"Verlet")

plt.scatter(dttab, x10)
plt.grid(True)
plt.xlabel("dt [s]")
plt.ylabel("x(10) [m]")
plt.show()

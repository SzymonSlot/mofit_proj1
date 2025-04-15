import numpy as np
import matplotlib.pyplot as plt
from funkcje import *
import csv

t = [0]
x = [2.8]
v = [0]
a = []
delta_t = []
e = 0.1
d=2
alfa = 0.5
tol = 1e-7
m=1
dt = 0.001

# jawny euler
dt = 0.001
m=1
tol = 1e-7
alfas = [0.5,5]
wynik_t = []
wynik_x = []
wynik_v = []
wynik_E = []
wynik_dt = []

for alfa in alfas:

    t = [0]
    x = [2.8]
    v = [0]
    a = []
    delta_t = []
    e = 0
    d=1
    while t[-1] <= 100:
        a.append(an(x[-1], m))

        # krok 2dt
        xk2, vk2 = trapez(x[-1], v[-1], 2*dt, m, alfa)

        # 2 kroki z dt
        xpk1, vpk1 = trapez(x[-1], v[-1], dt, m, alfa)
        xpk2, vpk2 = trapez(xpk1, vpk1, dt, m, alfa)

        ex = np.abs((xpk2-xk2)/(2**d-1))
        ev = np.abs((vpk2-vk2)/(2**d-1))

        if ex>ev:
            e = ex
        else:
            e = ev

        if e<=tol:
            x.append(xpk2)
            v.append(vpk2)
            t.append(t[-1]+2*dt)
            delta_t.append(2*dt)

        dt = 0.9*dt*(tol/e)**(1/(d+1))

    v = np.array(v)
    x = np.array(x)

    E = m*v**2/2 + V(x)

    wynik_t.append(t)
    wynik_x.append(x)
    wynik_v.append(v)
    wynik_E.append(E)
    wynik_dt.append(delta_t)


plotting(wynik_t, wynik_x, wynik_v, wynik_E, alfas, "")

# plt.plot(wynik_t[0][:-1], wynik_dt[0], label="alfa=0.5")
plt.scatter(wynik_t[1][:-1], wynik_dt[1], label="alfa=5")
plt.xlabel("Czas t [s]")
plt.ylabel("Dobrany krok czasowy $\Delta$t [s]")
plt.grid(True)
plt.legend()
plt.show()

#zapis do pliku
with open('alfa05.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['t', 'dt']) 
    for i in range(len(wynik_t[0][:-1])): 
        writer.writerow([wynik_t[0][i], wynik_dt[0][i]])

with open('alfa5.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['t', 'dt']) 
    for i in range(len(wynik_t[1][:-1])): 
        writer.writerow([wynik_t[1][i], wynik_dt[1][i]])
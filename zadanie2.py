import numpy as np
import matplotlib.pyplot as plt
from funkcje import *

# jawny euler
dt = 0.001
m=1
t = [0]
x = [2.8]
v = [0]
a = []
e = 0
d=1
tol = 1e-7
krokczasowy = []

while t[-1] <= 100:
    a.append(an(x[-1], m=m))

    # krok 2dt
    xk2 = x[-1] + v[-1]*2*dt
    vk2 = v[-1] + a[-1]*2*dt # u_k+2

    # 2 kroki z dt
    xpk1 = x[-1] + v[-1]*dt 
    vpk1 = v[-1] + a[-1]*dt # u'_k+1
    xpk2 = xpk1 + vpk1*dt
    vpk2 = vpk1 + an(xpk1)*dt # u'_k+2

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
        krokczasowy.append(2*dt)

    dt = 0.9*dt*(tol/e)**(1/(d+1))



print(len(t), tol)

v = np.array(v)
x = np.array(x)

E = m*v**2/2 + V(x)

wynik_t = [t]
wynik_x = [x]
wynik_v = [v]
wynik_E = [E]
wynik_krok = [krokczasowy]

# Verlet
m=1
t = [0]
x = [2.8]
v = [0]
e = 0
d=2
dt = 0.001
krokczasowy = []

while t[-1] <= 100:
    # krok 2dt
    xk2, vk2 = verlet(x[-1], v[-1], 2*dt, m)

    # 2 kroki z dt
    xpk1, vpk1 = verlet(x[-1], v[-1], dt, m)
    xpk2, vpk2 = verlet(xpk1, vpk1, dt, m)

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
        krokczasowy.append(2*dt)

    dt = 0.9*dt*(tol/e)**(1/(d+1))

print(len(t), tol)


v = np.array(v)
x = np.array(x)

E = m*v**2/2 + V(x)

wynik_t.append(t)
wynik_x.append(x)
wynik_v.append(v)
wynik_E.append(E)
wynik_krok.append(krokczasowy)

# RK4
dt = 0.001
m=1
t = [0]
x = [2.8]
v = [0]
a = []
e = 0
d=4
krokczasowy = []

while t[-1] <= 100:
    a.append(an(x[-1], m))

    # krok 2dt
    xk2, vk2 = RK4(x[-1], v[-1], 2*dt, m=m)

    # 2 kroki z dt
    xpk1, vpk1 = RK4(x[-1], v[-1], dt, m=m)
    xpk2, vpk2 = RK4(xpk1, vpk1, dt, m=m)

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
        krokczasowy.append(2*dt)

    dt = 0.9*dt*(tol/e)**(1/(d+1))



print(len(t), tol)

v = np.array(v)
x = np.array(x)

E = m*v**2/2 + V(x)

wynik_t.append(t)
wynik_x.append(x)
wynik_v.append(v)
wynik_E.append(E)
wynik_krok.append(krokczasowy)

plotting(wynik_t, wynik_x, wynik_v, wynik_E, ["Jawny Euler", "Verlet", "RK4"], "")

plt.plot(wynik_t[0][:-1], wynik_krok[0], label="Jawny Euler")
plt.plot(wynik_t[1][:-1], wynik_krok[1], label="Verlet")
plt.plot(wynik_t[2][:-1], wynik_krok[2], label="RK4")
plt.xlabel("Czas t [s]")
plt.ylabel("Dobrany krok czasowy $\Delta$t [s]")
plt.grid(True)
plt.legend()
plt.show()
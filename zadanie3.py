import numpy as np
import matplotlib.pyplot as plt
from funkcje import *
import pandas as pd

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
wynik_tRK4 = []
wynik_xRK4 = []
wynik_vRK4 = []
wynik_ERK4 = []
wynik_dtRK4 = []

for alfa in alfas:

    t = [0]
    x = [2.8]
    v = [0]
    a = []
    delta_t = []
    e = 0
    d=1
    while t[-1] <= 100:
        a.append(an(x[-1], valfa=v[-1]*alfa, m=m))

        # krok 2dt
        xk2 = x[-1] + v[-1]*2*dt
        vk2 = v[-1] + a[-1]*2*dt # u_k+2

        # 2 kroki z dt
        xpk1 = x[-1] + v[-1]*dt 
        vpk1 = v[-1] + a[-1]*dt # u'_k+1
        xpk2 = xpk1 + vpk1*dt
        vpk2 = vpk1 + an(xpk1, vpk1*alfa, m)*dt # u'_k+2

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

    t = [0]
    x = [2.8]
    v = [0]
    a = []
    delta_t = []
    e = 0
    d=4

    while t[-1] <= 100:
        a.append(an(x[-1], m))

        # krok 2dt
        xk2, vk2 = RK4(x[-1], v[-1], 2*dt, alfa, m)

        # 2 kroki z dt
        xpk1, vpk1 = RK4(x[-1], v[-1], dt, alfa, m)
        xpk2, vpk2 = RK4(xpk1, vpk1, dt, alfa, m)

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

    wynik_tRK4.append(t)
    wynik_xRK4.append(x)
    wynik_vRK4.append(v)
    wynik_ERK4.append(E)
    wynik_dtRK4.append(delta_t)



# plotting(wynik_t, wynik_x, wynik_v, wynik_E, alfas, "Dodatkowe opory ruchu Jawny Euler")
# plotting(wynik_tRK4, wynik_xRK4, wynik_vRK4, wynik_ERK4, alfas, "Dodatkowe opory ruchu RK4")
# plotting([wynik_t[0], wynik_tRK4[0]], [wynik_x[0], wynik_xRK4[0]], [wynik_v[0], wynik_vRK4[0]], [wynik_E[0], wynik_ERK4[0]], ["Jawny Euler", "RK4"], "Alfa = 0.5")
# plotting([wynik_t[1], wynik_tRK4[1]], [wynik_x[1], wynik_xRK4[1]], [wynik_v[1], wynik_vRK4[1]], [wynik_E[1], wynik_ERK4[1]], ["Jawny Euler", "RK4"], "Alfa = 5")

# żeby uzyskać wykres korku czasowego dla wszystkich trzech metod należy wykonać 
# skrypt zadanie4.py w którym zachodzi całkowanie ruchu wzorem trapezów

try:
    alfa05 = pd.read_csv("alfa05.csv")
    plt.plot(alfa05['t'], alfa05['dt'], label="Metoda trapezów")
except:
    print("Brakuje pliku alfa05.csv dla metody trapezów.")



plt.plot(wynik_t[0][:-1], wynik_dt[0], label="Jawny Euler")
plt.plot(wynik_tRK4[0][:-1], wynik_dtRK4[0], label="RK4")
plt.xlabel("Czas t [s]")
plt.ylabel("Dobrany krok czasowy $\Delta$t [s]")
plt.grid(True)
plt.legend()
plt.yscale("log")
plt.show()


try:
    alfa05 = pd.read_csv("alfa5.csv")
    plt.plot(alfa05['t'], alfa05['dt'], label="Metoda trapezów")
except:
    print("Brakuje pliku alfa5.csv dla metody trapezów.")


plt.plot(wynik_t[1][:-1], wynik_dt[1], label="Jawny Euler")
plt.plot(wynik_tRK4[1][:-1], wynik_dtRK4[1], label="RK4")
plt.xlabel("Czas t [s]")
plt.ylabel("Dobrany krok czasowy $\Delta$t [s]")
plt.grid(True)
plt.legend()
plt.xlim(0, 100)
# plt.ylim(0,1.1)
plt.yscale("log")
plt.show()
import numba as nb
import numpy as np
import matplotlib.pyplot as plt

@nb.njit
def V(x):
    l1 = 1
    l2 = 1/np.sqrt(8)
    return -np.exp(-x**2/l1**2) - 8*np.exp(-(x-2)**2/l2**2)

@nb.njit
def an(x, valfa=0, m=1, dx=0.001):
    return -(V(x + dx) - V(x - dx)) / (2 * m * dx) - valfa

@nb.njit
def verlet(xn, vn, dt, m=1):
    a1 = an(xn, m=m)
    x1 = xn + vn*dt+a1*dt**2/2
    a2 = an(x1, m=m)
    v1 = vn + dt/2*(a1+a2)

    return [x1, v1]

@nb.njit
def RK4(xn, vn, dt, alfa=0, m=1):
    k11 = vn
    k21 = an(xn, valfa=vn*alfa, m=m)
    k12 = vn + dt/2*k21
    k22 = an(xn + dt/2*k11, valfa=(vn+dt/2*k21)*alfa, m=m)
    k13 = vn + dt/2*k22
    k23 = an(xn + dt/2*k12, valfa=(vn+dt/2*k22)*alfa, m=m)
    k14 = vn + dt*k23
    k24 = an(xn + dt*k13, valfa=(vn+dt/2*k23)*alfa, m=m)
    return [xn + dt/6*(k11 + 2*k12 + 2*k13 + k14), 
            vn + dt/6*(k21 + 2*k22 + 2*k23 + k24)]


@nb.njit
def F1(xn, xn1, vn, vn1, dt):
    return xn1-xn-dt/2*vn1-dt/2*vn

@nb.njit
def F2(xn, xn1, vn, vn1, dt, alfa=0, m=1):
    return vn1-vn-dt/2*an(xn1, vn1*alfa, m=m)-dt/2*an(xn, vn*alfa, m=m)

@nb.njit
def dphi2dx2(xn, dx=0.001):
    return (V(xn+dx)-2*V(xn)+V(xn-dx))/dx**2


def trapez(xn, vn, dt, m=1, alfa=0):
    xmu = xn
    vmu = vn
    tol = 1e-7
    dxmu = 1
    dvmu = 1

    while(dxmu > tol and dvmu > tol):
        wsp = np.array([[1, -dt/2], [dt/(2*m)*dphi2dx2(xmu), 1+dt/2*alfa]])
        wynik = np.array([-F1(xn, xmu, vn, vmu, dt), -F2(xn, xmu, vn, vmu, dt, alfa, m)])

        # np.linalg.solve() rozwiązuje układ równań liwowych
        # dla układu x0 + 2 * x1 = 1 and 3 * x0 + 5 * x1 = 2
        # arr1 = np.array([[1, 2], [3, 5]])
        # arr2 = np.array([1, 2])
        # [x0, x1] = np.linalg.solve(arr1, arr2)

        dxmu, dvmu = np.linalg.solve(wsp, wynik)
        xmu = xmu + dxmu
        vmu = vmu + dvmu

    return [xmu, vmu]

def plotting(t_list, x_list, v_list, E_list, dt_list, title=""):
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    colors = ["skyblue", "fuchsia", "red", "orange", "purple", "lime"]
    
    for i, (t, x) in enumerate(zip(t_list, x_list)):
        axs[0, 0].plot(t, x, label=f"{dt_list[i]}", color=colors[i % len(colors)])
    axs[0, 0].set_xlabel("t [s]")
    axs[0, 0].set_ylabel("x [m]")
    axs[0, 0].set_title("Położenie w czasie")
    axs[0, 0].legend()
    axs[0, 0].grid()
    # axs[0, 0].set_xlim(0, 20)

    for i, (t, v) in enumerate(zip(t_list, v_list)):
        axs[0, 1].plot(t, v, label=f"{dt_list[i]}", color=colors[i % len(colors)])
    axs[0, 1].set_xlabel("t [s]")
    axs[0, 1].set_ylabel("v [m/s]")
    axs[0, 1].set_title("Prędkość w czasie")
    axs[0, 1].legend()
    axs[0, 1].grid()
    # axs[0, 1].set_xlim(0, 20)
    
    for i, (x, v) in enumerate(zip(x_list, v_list)):
        axs[1, 0].plot(x, v, label=f"{dt_list[i]}", color=colors[i % len(colors)])
    axs[1, 0].set_xlabel("x [m]")
    axs[1, 0].set_ylabel("v [m/s]")
    axs[1, 0].set_title("Portret fazowy")
    axs[1, 0].legend()
    axs[1, 0].grid()

    for i, (t, E) in enumerate(zip(t_list, E_list)):
        axs[1, 1].plot(t, E, label=f"{dt_list[i]}", color=colors[i % len(colors)])
    axs[1, 1].set_xlabel("t [s]")
    axs[1, 1].set_ylabel("E [J]")
    axs[1, 1].set_title("Energia w czasie")
    axs[1, 1].legend()
    axs[1, 1].grid()
    # axs[1, 1].set_xlim(0, 20)
    
    plt.tight_layout()
    plt.show()
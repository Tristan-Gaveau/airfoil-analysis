import numpy as np
import matplotlib.pyplot as plt
from geometry import naca4_sym
from panels import resoudre_kutta

# XFOIL : alpha en colonne 0, CL en colonne 1
d = np.loadtxt('results/xfoil_polaire.txt', skiprows=12)
alpha_x, Cl_x = d[:, 0], d[:, 1]

# Ton code, aux mêmes angles
x, yt = naca4_sym(n=100)
xf = np.concatenate((x[::-1], x[1:]))
yf = np.concatenate((-yt[::-1], yt[1:]))

Cl_mon = []
for a in alpha_x:
    al = np.radians(a)
    lam, gamma, Vt, Cp, xc, yc, ds, nx, ny = resoudre_kutta(xf, yf, al)
    Fx = -np.sum(Cp * nx * ds)
    Fy = -np.sum(Cp * ny * ds)
    Cl_mon.append(Fy*np.cos(al) - Fx*np.sin(al))
Cl_mon = np.array(Cl_mon)

print(" alpha   Cl_XFOIL   Cl_mon   ecart(abs)")
for a, c1, c2 in zip(alpha_x, Cl_x, Cl_mon):
    print(f"{a:6.1f}  {c1:9.4f}  {c2:8.4f}  {c2 - c1:+9.4f}")

plt.plot(alpha_x, Cl_x, 's-', label='XFOIL (non visqueux)')
plt.plot(alpha_x, Cl_mon, 'o--', label='Mon code')
plt.xlabel('alpha (°)'); plt.ylabel('Cl'); plt.grid(True); plt.legend()
plt.savefig('results/comparaison_polaire.png', dpi=150)
plt.show()
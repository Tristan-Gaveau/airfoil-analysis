import numpy as np
import matplotlib.pyplot as plt
from geometry import naca4_sym
from panels import resoudre_kutta

x, yt = naca4_sym(n=100)
xf = np.concatenate((x[::-1], x[1:]))
yf = np.concatenate((-yt[::-1], yt[1:]))

alphas = np.arange(-4, 9, 2)
Cls = []
for a in alphas:
    al = np.radians(a)
    lam, gamma, Vt, Cp, xc, yc, ds, nx, ny = resoudre_kutta(xf, yf, al)
    Fx = -np.sum(Cp * nx * ds)
    Fy = -np.sum(Cp * ny * ds)
    Cl = Fy*np.cos(al) - Fx*np.sin(al)
    Cls.append(Cl)
    print(f"alpha={a:3d}  Cl={Cl:7.4f}  KJ={-2*gamma*ds.sum():7.4f}  Kutta={Vt[0]+Vt[-1]:.1e}")

pente = np.polyfit(np.radians(alphas), Cls, 1)[0]
print("Pente dCl/dalpha :", pente, "(2*pi =", 2*np.pi, ")")

plt.plot(alphas, Cls, 'o-')
plt.xlabel('alpha (°)'); plt.ylabel('Cl'); plt.grid(True)
plt.savefig('results/naca0012_cl_alpha.png', dpi=150)
plt.show()
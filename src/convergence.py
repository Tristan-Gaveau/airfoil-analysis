import numpy as np
import matplotlib.pyplot as plt
from geometry import naca4_sym
from panels import resoudre_kutta

alpha = np.radians(5.0)
valeurs_n = [25, 50, 100, 200]
nb_panneaux = []
Cls = []

for n in valeurs_n:
    # Géométrie
    x, yt = naca4_sym(n=n)
    xf = np.concatenate((x[::-1], x[1:]))
    yf = np.concatenate((-yt[::-1], yt[1:]))

    # Calcul
    lam, gamma, Vt, Cp, xc, yc, ds, nx, ny = resoudre_kutta(xf, yf, alpha)

    # Forces
    Fx = -np.sum(Cp * nx * ds)
    Fy = -np.sum(Cp * ny * ds)
    Cl = Fy*np.cos(alpha) - Fx*np.sin(alpha)

    nb_panneaux.append(len(xc))
    Cls.append(Cl)
    print(f"n={n:4d}  panneaux={len(xc):4d}  Cl={Cl:.4f}")

# Écart par rapport au cas le plus fin
Cl_ref = Cls[-1]
for npan, Cl in zip(nb_panneaux, Cls):
    print(f"{npan:4d} panneaux : écart = {100*(Cl - Cl_ref)/Cl_ref:+.3f} %")

plt.plot(nb_panneaux, Cls, 'o-')
plt.xlabel('Nombre de panneaux')
plt.ylabel('Cl (alpha = 5°)')
plt.title('Convergence de Cl en fonction du nombre de panneaux')
plt.grid(True)
plt.savefig('results/convergence_cl.png', dpi=150)
plt.show()
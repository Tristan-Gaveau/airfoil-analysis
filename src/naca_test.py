import numpy as np
import matplotlib.pyplot as plt
from geometry import naca4_sym
from panels import resoudre_sources, contours_panneaux

alpha_deg = 0.0                  # change cette valeur : 0.0, puis 5.0
alpha = np.radians(alpha_deg)

# 1. Géométrie
x, yt = naca4_sym(n=100)
xf = np.concatenate((x[::-1], x[1:]))
yf = np.concatenate((-yt[::-1], yt[1:]))   # sens horaire, bord de fuite fermé

# 2. Calcul
lam, Vt, Cp, xc, yc, ds = resoudre_sources(xf, yf, alpha)
_, _, _, _, _, _, _, nx, ny = contours_panneaux(xf, yf)

# 3. Forces
Fx = -np.sum(Cp * nx * ds)
Fy = -np.sum(Cp * ny * ds)
Cl = Fy*np.cos(alpha) - Fx*np.sin(alpha)
Cd = Fy*np.sin(alpha) + Fx*np.cos(alpha)

# 4. Affichage
print("Nombre de panneaux :", len(xc))
print("Flux total :", np.sum(lam*ds))
print("Cp min :", Cp.min(), "à x =", xc[Cp.argmin()])
print("Cp max :", Cp.max())
print("Cl =", Cl, " Cd =", Cd)

plt.plot(xc, Cp, '.-')
plt.gca().invert_yaxis()          # convention aérodynamique : Cp négatif en haut
plt.xlabel('x/c'); plt.ylabel('Cp'); plt.grid(True)
plt.title(f'NACA 0012, alpha = {alpha_deg}°, sources seules')
plt.savefig(f'results/naca0012_alpha{alpha_deg:g}.png', dpi=150)
plt.show()
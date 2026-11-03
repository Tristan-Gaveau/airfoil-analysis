import numpy as np

# Air standard au niveau de la mer
T = 288.15        # K
rho = 1.225       # kg/m3
mu = 1.7894e-5    # Pa.s
gamma, R = 1.4, 287.05

M = 0.15
Re = 2e6
y_plus = 1.0

a = 340.3        # vitesse du son
V = M * a         # vitesse
c = 0.5723        # corde
q = 1/2 * rho * V**2        # pression dynamique
Cf = 0.00327       # coefficient de frottement estimé
tau_w = Cf * q    # contrainte à la paroi; formule : Cf * q
u_tau = np.sqrt(tau_w / rho)    # vitesse de frottement; formule : np.sqrt(tau_w / rho)
y1 = y_plus * mu / (rho * u_tau)      # hauteur de la première cellule
delta = 11.6e-3    # épaisseur de couche limite

print(f"Vitesse du son : {a:.1f} m/s")
print(f"Vitesse V      : {V:.2f} m/s")
print(f"Corde c        : {c:.4f} m")
print(f"Pression dyn.  : {q:.1f} Pa")
print(f"Cf estimé      : {Cf:.5f}")
print(f"u_tau          : {u_tau:.3f} m/s")
print(f"y1 (y+=1)      : {y1*1e6:.2f} micromètres")
print(f"delta          : {delta*1000:.2f} mm")
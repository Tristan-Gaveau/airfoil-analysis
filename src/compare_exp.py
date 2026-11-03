import numpy as np
import matplotlib.pyplot as plt
from geometry import naca4_sym
from panels import resoudre_kutta

# --- 1. Ton code, de 0° à 20°
x, yt = naca4_sym(n=100)
xf = np.concatenate((x[::-1], x[1:]))
yf = np.concatenate((-yt[::-1], yt[1:]))

alphas = np.arange(0, 21, 2)
Cl_mon = []
for a in alphas:
    al = np.radians(a)
    lam, gamma, Vt, Cp, xc, yc, ds, nx, ny = resoudre_kutta(xf, yf, al)
    Fx = -np.sum(Cp * nx * ds)
    Fy = -np.sum(Cp * ny * ds)
    Cl_mon.append(Fy*np.cos(al) - Fx*np.sin(al))
Cl_mon = np.array(Cl_mon)

plt.plot(alphas, Cl_mon, 'o--', label='Mon code (non visqueux)')

# --- 2. Les deux courbes expérimentales
fichiers = {
    'Ladson Re = 2e6': 'results/exp_Re2e6.csv',
    'Ladson Re = 11.9e6': 'results/exp_Re12e6.csv',
}

for nom, chemin in fichiers.items():
    exp = np.loadtxt(chemin, delimiter=',', skiprows=1)
    alpha_e, Cl_e = exp[:, 0], exp[:, 1]
    Cl_mon_aux_exp = np.interp(alpha_e, alphas, Cl_mon)

    print(nom)
    print(" alpha   Cl_exp   Cl_mon   ecart(abs)   ecart(%)")
    for a, ce, cm in zip(alpha_e, Cl_e, Cl_mon_aux_exp):
        pct = 100*(cm - ce)/ce if ce != 0 else float('nan')
        print(f"{a:6.1f}  {ce:7.3f}  {cm:7.3f}  {cm - ce:+9.3f}  {pct:+8.1f}")
    plt.plot(alpha_e, Cl_e, 's-', label=nom)

# --- 3. Graphique
plt.xlabel('alpha (°)')
plt.ylabel('Cl')
plt.grid(True)
plt.legend()
plt.title('NACA 0012 : code vs expérience')
plt.savefig('results/comparaison_exp.png', dpi=150)
plt.show()
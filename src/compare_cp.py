import numpy as np
import matplotlib.pyplot as plt
from geometry import naca4_sym
from panels import resoudre_kutta

alpha = np.radians(5.0)

# --- XFOIL : colonnes x, y, Cp. Le fichier part du bord de fuite, passe par
# --- l'extrados jusqu'au bord d'attaque, puis revient par l'intrados.
d = np.loadtxt('results/cp_xfoil_5deg.txt', skiprows=3)
x_x, cp_x = d[:, 0], d[:, 2]
iLE = np.argmin(x_x)                              # bord d'attaque = x minimal
xu, cu = x_x[:iLE+1][::-1], cp_x[:iLE+1][::-1]    # extrados, x croissant
xl, cl = x_x[iLE:], cp_x[iLE:]                    # intrados, x croissant

# --- Ton calcul : sens horaire, donc d'abord l'intrados (bord de fuite -> bord d'attaque)
x, yt = naca4_sym(n=100)
xf = np.concatenate((x[::-1], x[1:]))
yf = np.concatenate((-yt[::-1], yt[1:]))
lam, gamma, Vt, Cp, xc, yc, ds, nx, ny = resoudre_kutta(xf, yf, alpha)

demi = len(xc) // 2
xc_l, Cp_l = xc[:demi], Cp[:demi]                 # intrados
xc_u, Cp_u = xc[demi:], Cp[demi:]                 # extrados

# --- XFOIL interpolé sur tes points de contrôle
ref_u = np.interp(xc_u, xu, cu)
ref_l = np.interp(xc_l, xl, cl)

# --- Écarts, hors bords d'attaque et de fuite (le Cp y varie trop vite)
m_u = (xc_u > 0.02) & (xc_u < 0.98)
m_l = (xc_l > 0.02) & (xc_l < 0.98)
ecart_u = np.abs(Cp_u - ref_u)[m_u]
ecart_l = np.abs(Cp_l - ref_l)[m_l]
print("Extrados : écart moyen", ecart_u.mean(), " max", ecart_u.max())
print("Intrados : écart moyen", ecart_l.mean(), " max", ecart_l.max())

plt.plot(x_x, cp_x, '-', label='XFOIL')
plt.plot(xc, Cp, '.', label='Mon code')
plt.gca().invert_yaxis()
plt.xlabel('x/c'); plt.ylabel('Cp'); plt.grid(True); plt.legend()
plt.savefig('results/comparaison_cp_5deg.png', dpi=150)
plt.show()
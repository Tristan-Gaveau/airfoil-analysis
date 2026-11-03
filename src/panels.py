import numpy as np
import matplotlib.pyplot as plt

def contours_panneaux(x, y):
    x1, y1 = x[:-1], y[:-1]
    x2, y2 = x[1:], y[1:]

    dx = x2 - x1
    dy = y2 - y1
    xc = (x1 + x2) / 2
    yc = (y1 + y2) / 2
    ds = np.sqrt(dx**2 + dy**2)

    nx = -dy / ds    # normale extérieure (sens horaire)
    ny = dx / ds

    return x1, y1, x2, y2, ds, xc, yc, nx, ny

def influence(xp, yp, x1, y1, x2, y2):
    """Vitesse (u, v) induite en (xp, yp) par un panneau de source unitaire."""

    # Géométrie du panneau
    dx = x2 - x1
    dy = y2 - y1
    ds = np.sqrt(dx**2 + dy**2)
    phi = np.arctan2(dy, dx)

    # Étape 1 : repère global -> repère du panneau
    xi = (xp - x1)*np.cos(phi) + (yp - y1)*np.sin(phi)
    eta = -(xp - x1)*np.sin(phi) + (yp - y1)*np.cos(phi)

    # Étape 2 : vitesse dans le repère du panneau
    r1_sq = xi**2 + eta**2
    r2_sq = (xi - ds)**2 + eta**2
    theta1 = np.arctan2(eta, xi)
    theta2 = np.arctan2(eta, xi - ds)
    u_xi = (1/(4*np.pi)) * np.log(r1_sq / r2_sq)
    u_eta = (1/(2*np.pi)) * (theta2 - theta1)

    # Étape 3 : retour au repère global
    u = u_xi*np.cos(phi) - u_eta*np.sin(phi)
    v = u_xi*np.sin(phi) + u_eta*np.cos(phi)

    return u, v

def resoudre_sources(x, y, alpha=0.0):
    x1, y1, x2, y2, ds, xc, yc, nx, ny = contours_panneaux(x, y)
    N = len(xc)
    tx, ty = ny, -nx

    A = np.zeros((N, N))
    T = np.zeros((N, N))      # influence tangentielle
    for i in range(N):
        for j in range(N):
            if i == j:
                A[i, j] = 0.5
                T[i, j] = 0.0
            else:
                u, v = influence(xc[i], yc[i], x1[j], y1[j], x2[j], y2[j])
                A[i, j] = u * nx[i] + v * ny[i]
                T[i, j] = u * tx[i] + v * ty[i]

    b = -(np.cos(alpha) * nx + np.sin(alpha) * ny)
    lam = np.linalg.solve(A, b)

    Vt = np.cos(alpha) * tx + np.sin(alpha) * ty + T @ lam
    Cp = 1 - Vt**2
    return lam, Vt, Cp, xc, yc, ds

def resoudre_kutta(x, y, alpha=0.0):
    x1, y1, x2, y2, ds, xc, yc, nx, ny = contours_panneaux(x, y)
    N = len(xc)
    tx, ty = ny, -nx

    An = np.zeros((N, N)); At = np.zeros((N, N))   # sources
    Bn = np.zeros((N, N)); Bt = np.zeros((N, N))   # vortex
    for i in range(N):
        for j in range(N):
            if i == j:
                An[i, j] = 0.5;  At[i, j] = 0.0
                Bn[i, j] = 0.0;  Bt[i, j] = -0.5
            else:
                u, v = influence(xc[i], yc[i], x1[j], y1[j], x2[j], y2[j])
                An[i, j] = u*nx[i] + v*ny[i]     # composante normale
                At[i, j] = u*tx[i] + v*ty[i]     # composante tangentielle
                uv, vv = -v, u
                Bn[i, j] = uv*nx[i] + vv*ny[i]   # composante normale du vortex
                Bt[i, j] = uv*tx[i] + vv*ty[i]   # composante tangentielle du vortex

    M = np.zeros((N+1, N+1))
    b = np.zeros(N+1)
    M[:N, :N] = An
    M[:N, N] = Bn.sum(axis=1)
    b[:N] = -(np.cos(alpha) * nx + np.sin(alpha) * ny)
    M[N, :N] = At[0] + At[N-1]
    M[N, N] = Bt[0].sum() + Bt[N-1].sum()
    b[N] = -(np.cos(alpha)*(tx[0]+tx[N-1]) + np.sin(alpha)*(ty[0]+ty[N-1]))

    sol = np.linalg.solve(M, b)
    lam, gamma = sol[:N], sol[N]
    Vt = np.cos(alpha)*tx + np.sin(alpha)*ty + At @ lam + gamma*Bt.sum(axis=1)    
    Cp = 1 - Vt**2
    return lam, gamma, Vt, Cp, xc, yc, ds, nx, ny

if __name__ == "__main__":
    N = 100
    theta_pts = np.linspace(2*np.pi, 0, N+1)   # sens horaire
    x, y = np.cos(theta_pts), np.sin(theta_pts)

    # --- Résolution ---
    lam, Vt, Cp, xc, yc, ds = resoudre_sources(x, y)

    # --- Solutions exactes ---
    theta = np.arctan2(yc, xc)
    Cp_exact = 1 - 4*np.sin(theta)**2
    lam_exact = -2*np.cos(theta)

    # --- Les trois tests ---
    print("Flux total (attendu ~ 0) :", np.sum(lam*ds))
    print("Écart max Cp :", np.max(np.abs(Cp - Cp_exact)))
    print("Écart max lambda :", np.max(np.abs(lam - lam_exact)))

    # --- Tracé ---
    ordre = np.argsort(theta)
    plt.plot(np.degrees(theta[ordre]), Cp[ordre], 'o', label='Panneaux')
    plt.plot(np.degrees(theta[ordre]), Cp_exact[ordre], '-', label='Exact')
    plt.xlabel('θ (°)')
    plt.ylabel('Cp')
    plt.legend()
    plt.grid(True)
    plt.savefig('results/cylindre_cp.png', dpi=150)
    plt.show()
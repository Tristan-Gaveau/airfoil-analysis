import numpy as np
from scipy.integrate import quad
from panels import contours_panneaux

def systeme_sources(x, y, alpha_deg):
    x1, y1, x2, y2, ds, xc, yc, nx, ny = contours_panneaux(x, y)
    N = len(ds)
    alpha = np.radians(alpha_deg)
    A = np.zeros((N, N))
    b = np.zeros(N)

    for i in range(N):
        b[i] = -(np.cos(alpha) * nx[i] + np.sin(alpha) * ny[i])
        for j in range(N):
            if i == j:
                A[i, j] = 0.5
            else:
                def integrande(t):
                    xi = x1[j] + t * (x2[j] - x1[j])
                    eta = y1[j] + t * (y2[j] - y1[j])
                    num = ((xc[i] - xi) * nx[i] + (yc[i] - eta) * ny[i])           # (xc_i - xi)·nx_i + (yc_i - eta)·ny_i
                    den = ((xc[i] - xi)**2 + (yc[i] - eta)**2)           # (xc_i - xi)² + (yc_i - eta)²
                    return num / den * ds[j]
                I, _ = quad(integrande, 0, 1)
                A[i, j] = I / (2 * np.pi)
    return A, b, ds, nx, ny


if __name__ == "__main__":
    N = 50
    theta = np.linspace(2*np.pi, 0, N + 1)
    x, y = np.cos(theta), np.sin(theta)

    A, b, ds, nx, ny = systeme_sources(x, y, alpha_deg=5)
    lam = np.linalg.solve(A, b)

    print("Somme lam * ds :", np.sum(lam * ds))
    print("Écart max avec 2*b :", np.max(np.abs(lam - 2 * b)))
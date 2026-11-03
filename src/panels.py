import numpy as np

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


if __name__ == "__main__":
    N = 50
    theta = np.linspace(2*np.pi, 0, N+1)   # sens horaire
    x, y = np.cos(theta), np.sin(theta)

    x1, y1, x2, y2, ds, xc, yc, nx, ny = contours_panneaux(x, y)

    print("Somme des ds :", ds.sum(), "(attendu ~ 6.279)")
    print("Écart max normale :", np.max(np.abs(nx - xc)), np.max(np.abs(ny - yc)))
import numpy as np
import matplotlib.pyplot as plt

def naca4_sym(t=0.12, n=100):
    beta = np.linspace(0, np.pi, n) # n angles entre 0 et pi
    x = 0.5 * (1 - np.cos(beta)) # points de la corde de l'aile
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4) # la demi épaisseur de l'aile (formule de NACA)
    return x, yt #renvoie les positions x et la demi-épaisseurs yt

if __name__ == "__main__":
    x, yt = naca4_sym()
    plt.plot(x, yt, 'b')
    plt.plot(x, -yt, 'b')
    plt.axis('equal')
    plt.show()
    print("Épaisseur max :", 2 * yt.max())
    print("Position :", x[yt.argmax()])
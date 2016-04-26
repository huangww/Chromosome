import numpy as np
import matplotlib.pyplot as plt

def flj(d0):
    r0 = 0.5
    eps = 10.0
    if d0**6 <= 2*r0**6:
        d6 = (r0/d0)**6
        f = 48*eps*(d6-0.5)*d6/d0
    else:
        f = 0
    return f

def main():
    x = np.linspace(0.001, 1, 100)
    f = np.array([flj(d) for d in x])
    plt.plot(x, f)
    plt.ylim([-1, 100])
    plt.show()

if __name__ == "__main__":
    main()



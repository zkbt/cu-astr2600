
import numpy as np

default = 50
a = 10
b = -10.0
c = 1.0
sigma = 4.0
def flat(N=default, seed=42):
    np.random.seed(seed)
    x = np.linspace(-5, 5, N)
    m = a
    y = m + np.random.normal(0, sigma, size=N)
    return x, y, sigma

def line(N=default, seed=42):
    np.random.seed(seed)
    x = np.linspace(-5, 5, N)
    m = a + b*x
    y = m + np.random.normal(0, sigma, size=N)
    return x, y, sigma

def parabola(N=default, seed=42):
    np.random.seed(seed)
    x = np.linspace(-5, 5, N)
    m = a + b*x + c*x**2
    y = m + np.random.normal(0, sigma, size=N)
    return x, y, sigma

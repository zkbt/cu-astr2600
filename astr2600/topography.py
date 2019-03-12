import numpy as np

def simulate_topography(numCraters = 10, # number of craters
                        widthMax   = 1.0,     # maximal width of Gaussian crater
                        heightMin  = -1.0,    # maximal depth of craters / valleys
                        heightMax  = 2.0):     # maximal height of hills / mountains

    # 1-D Gaussian
    np.random.seed()   # for reproducability
    def gaussian(x, A, mu, sigma):
        return A * np.exp(-(x - mu)**2 / 2.0 / sigma**2)

    # 1-D Gaussian (same thing using lambda)
    #gaussian = lambda x, A, mu, sigma: A * np.exp(-(x - mu)**2 / 2. / sigma**2)

    # Create an array of linearly spaced x values
    xArr = np.linspace(0, 10, 500)  # km

    # Create an array of initially flat landscape (aka filled with 0's)
    yArr = np.zeros_like(xArr)

    # Add craters / mountains to landscape
    for _ in range(numCraters):

        # Amplitude between heightMin and heightMax
        A = np.random.rand() * (heightMax - heightMin) + heightMin

        # Center location of the crater
        center = np.random.rand() * xArr.max()

        # Width of the crater
        sigma = np.random.rand() * widthMax

        # Add crater to landscape!
        yArr += gaussian(xArr, A=A, mu=center, sigma=sigma)

    return xArr, yArr

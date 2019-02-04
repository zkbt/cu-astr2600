import numpy as np
import matplotlib.pyplot as plt

def drawMoons(names, xpositions, xlim=[-500,500], labels=False):
    '''
    Draw a plot of the positions of moons relative to Jupiter.

    This function requires two input arguments. They should
    both be lists, and they should be the same size as each
    other. They are:

        moons -- a 1-dimensional list of moon names
        xpositions -- a 1-dimensional list of moon positions (in arcsec)

    For example:

        %matplotlib inline
        names = ['Io', 'Europa', 'Ganymede', 'Callisto']
        xpositions = [-20, 40, 80, -160]
        drawMoons(names, xpositions)

        (should display a plot of the moon positions)

    Options keyword arguments

        xlim = [-500,500]
            (set this to some other 2-element list to change
                the range over which the plot stretchs)
        labels = False
            (set this to True to provide text labels)
    '''

    # since we're plotting only 1D positions, we make up y-values
    ypositions = np.zeros_like(xpositions)

    # we create a new figure, and set its size
    plt.figure(figsize=(10,.25))

    # we plot the moons in their positions
    plt.plot(xpositions, ypositions,
                marker = '.',
                linewidth=0,
                color='black')

    # if desired, we add text labels to all the moons
    if labels:
        for x, y, n in zip(xpositions, ypositions, names):
            plt.text(x, y+0.1, n, ha='center', va='bottom')

    # plot Jupiter in the center
    plt.plot(0,0, marker='o', markersize=10, color='black')

    # set the x and y limits of the plot
    plt.xlim(*xlim)
    plt.ylim(-1,1)

    # turn off all axis labels (and the box around the plot)
    plt.axis('off')

    # make sure the plot shows to the screen
    plt.show()

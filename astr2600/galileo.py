import numpy as np
import matplotlib.pyplot as plt

def drawMoons(names, xpositions, xlim=[-500,500], labels=True):
    '''
    Draw a plot of the positions of moons relative to Jupiter.

    This function requires two input arguments. They should
    both be lists, and they should be the same size as each
    other. They are:

        moons -- a 1-dimensional list of moon names
        xpositions -- a 1-dimensional list of moon positions (in arcsec)

    For example:

        names = ['Io', 'Europa', 'Ganymede', 'Callisto']
        xpositions = [-20, 40, 80, -160]
        drawMoons(names, xpositions)

        (this should display a plot of the moon positions)

    Options keyword arguments

        xlim = [-500,500]
            This defines the x values of the left and
            right edges of the plotting range to be included.

        labels = True
            If the function is called with labels=True,
            then display the names of the moons.

            If the function is called with labels=False,
            then do not display the names of the moons.
    '''

    # since we're plotting only 1D positions, we make up y-values
    ypositions = np.zeros_like(xpositions)

    # we create a new figure, and set its size
    plt.figure(figsize=(10,0.5))

    # we plot the moons in their positions
    plt.plot(xpositions, ypositions,
                marker = '.',
                linewidth=0,
                color='black')

    # if desired, we add text labels to all the moons
    if labels:
        for x, y, n in zip(xpositions, ypositions, names):
            plt.text(x, y+0.5, n, ha='center', va='bottom', size=9)

    # plot Jupiter in the center
    plt.plot(0,0, marker='o', markersize=20, markerfacecolor='none', markeredgecolor='black')

    # set the x and y limits of the plot
    plt.xlim(*xlim)
    plt.ylim(-1,1)

    # turn off all axis labels (and the box around the plot)
    plt.axis('off')

    # make sure the plot shows to the screen
    plt.show()

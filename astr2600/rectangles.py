'''
This module contains tools for estimating numerical integrals.
'''

import numpy as np
import matplotlib.pyplot as plt

def visualizeRectangles(f, xmin, xmax, N):
    """

    Plot "N" rectangles between "xmin" and "xmax",
    to visualize what happens when we calculate
    a rectangular integral of the function "f".

    Arguments
        f    : a function to be integrated
        xmin : lower limit of integration
        xmax : upper limit of integration
        N    : number of rectangles to use

    Returns
        (nothing), but it makes a plot

    Example:
        def integrand(x)
            return x**2
        plotRectangles(integrand, 0, 2, 5)
    """

    # first, make a high-resolution plot of the function
    xsmooth = np.linspace(xmin, xmax, 500)
    plt.plot(xsmooth, f(xsmooth), linewidth=10, alpha=0.3, color='orchid')

    # next, make a grid spanning xmin to xmax, with N+1 elements
    x = np.linspace(xmin, xmax, N+1)

    # calculate the offsets from x[:-1] to x[1:]
    dx = np.diff(x)

    # loop over rectangles
    for i in range(N):

        # where are the left and right edges of the rectangle?
        left = x[i]
        right = x[i] + dx[i]

        # set the rectangle plotting style
        plotkw = dict(color='black', linewidth=2)

        # plot vertical line from 0 up to curve (left edge of rectangle)
        plt.plot([left, left], [0, f(x[i])], **plotkw)

        # plot another vertical line (right edge of rectangle)
        plt.plot([right, right], [0, f(x[i])], **plotkw)

        # plot horizontal line (top of rectangle)
        plt.plot([left, right], [f(x[i]), f(x[i])], **plotkw)

    # set the x and y labels
    plt.xlabel('x')
    plt.ylabel('f(x)')

    # make sure the plot displays to the screen
    plt.show()

def integrateRectangles(f, xmin, xmax, N):
    """

    Approximate the area under the function "f",
    from "xmin" to "xmax", using N rectangles.

    Arguments
        f    : a function to be integrated
        xmin : lower limit of integration
        xmax : upper limit of integration
        N    : number of rectangles to use

    Returns
        the numerical integral of "f"

    Example:
        def integrand(x)
            return x**2
        integrateRectangles(integrand, 0, 1, 5)
    """

    # These x values represent edges of rectangles.
    # We need (N+1) edges for N rectangles.
    x = np.linspace(xmin, xmax, N + 1)

    # The widths of the rectangles:
    dx = np.diff(x)  # equivalently x[1:] - x[:-1]

    # The heights of the rectangles (using the left of each pair of points)
    h = f(x[:-1])

    # The areas of each rectangle.
    A = dx * h

    # Add up all those areas.
    return np.sum(A)

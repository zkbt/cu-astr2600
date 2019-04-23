import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

h = 6.626068e-34    # Planck's constant, m**2 kg / s
k = 1.3806503e-23   # Boltzmann's constant, m**2 kg / s**2 / K
c = 299792458.0     # speed of light, m / s


class ThingThatGlows:
    '''
    The class of ThingThatGlows python objects can handle calculations
    and visualizations related to physical things that emit radiation.'''

    def __init__(self, radius, temperature):
        '''
        Initialize a ThingThatGlows object.

        Parameters
        ----------
        radius : float
            The radius of the object, in meters
        temperature : float
            The temperature of the object, in K
        '''

        # store the radius and temperature as attributes
        self.radius = radius
        self.temperature = temperature

    def __repr__(self):
        '''
        If we represent this object as a string,
        what will it look like?
        '''
        return '{} (R={}m, T={}K)'.format(self.__class__.__name__,
                                        self.radius,
                                        self.temperature)


    def calculateArea(self):
        '''
        Calculate the surface area of the thing,
        in units of meters**2.
        '''
        return 4 * np.pi * self.radius**2

    def calculateLuminosity(self):
        '''
        Calculate the luminosity of the thing,
        in units of W.
        '''

        # the Stefan-Bolztmann constant
        sigma = 5.670367e-08 # W/m**2/K**4

        # calculate the surface flux, in units of W/m**2
        # powerPerArea = sigma * self.temperature**4
        # integrate the flux from 100 to 10000
        def logf(logw):
            w = np.exp(logw)
            return w*self.flux(w)

        powerPerArea = quad(logf, np.log(100), np.log(100000))[0]

        # calculate the luminosity as flux * area
        return powerPerArea * self.calculateArea()

    def plotSpectrum(self):
        '''
        Plot the spectrum of light this thing emits.
        '''

        # create a grid of wavelengths
        wavelengths = np.logspace(2, 6, 1000)

        # this is the flux per unit area
        flux = self.flux(wavelengths)

        # we multiply by an area to get total luminosity
        specificluminosity = flux * self.calculateArea()

        # make the plot
        plt.plot(wavelengths, specificluminosity, label=repr(self))
        plt.xscale('log'); plt.yscale('log')
        plt.xlabel('Wavelength (nm)'); plt.ylabel('Specific Luminosity (W/nm)')
        plt.legend(frameon=False)

    def flux(self, w_nm):
        '''
        Compute the thermal emission spectrum.

        Parameters
        ----------
        w_nm : np.array
            Wavelengths, in units of nanometer


        Returns
        -------
        Flux at the given wavelength(s), in W / m**2 / nm
        '''

        # convert nm to m
        w = w_nm/1e9

        # pull out the temperature from this object's attribute
        T = self.temperature

        # return Planck function, in W/m**2/nm/sr
        num = 2 * h *c**2
        den = w**5 * (np.exp(h * c / (w * k * T)) - 1.0)
        return np.pi * num / den / 1e9

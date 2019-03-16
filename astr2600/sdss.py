'''
This module contains tools to download astronomical spectra from SDSS.

Functions:

    getSDSS: Donwload, extract, and return a Sloan Digital Sky Survey (SDSS) spectrum from the web

    getLabel: Returns a metadata label (string) corresponding to a given spectrum obtained from getSDSS

    test: A function to test the connection to the SDSS Server, and either download a
        spectrum or generate a fake spectrum

    fakeSDSS: Generate and return a fake SDSS spectrum (a random-temperature blackbody spectrum)

    planck:  Calculate and return the Planck function for a given temperature

Constants:
    All constants given in SI units

    h: Planck's constant

    k: Boltzmann's Constant

    c: Speed of Light
'''

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pf


def downloadSDSS(online=True, verbose=True):
    '''
    Return a random spectrum from the Sloan Digital Sky Survey (SDSS).

        Arguments:

            online:  Boolean keyword to enable testing offline *if* the
                spectrum spec-4055-55359-0596.fits is already downloaded
            verbose: Boolean keyword:  if False, omit print statements

        Returns:

            one random SDSS spectrum. The format of this spectrum
            is as a dictionary that contains many key-value pairs.
            Important keys inside this spectrum are highlighted in
            the example below.

        Example:

            # call the function to get a random spectrum
            s = downloadSDSS()

            # the wavelengths are in units of angstroms (=10**-10 m)
            w_angstroms = s['wavelengths']
            w_nm = np.array(w_angstroms)/10.0

            # the fluxes are in units of "1e-17 erg/cm**2/s/angstrom"
            f = s['flux']

    '''

    # choose spectrum to download
    plate = '4055'
    mjd = '55359'
    num = np.random.randint(0,1000)  # choose among the 1000 spectra here
    numstr = str(num)
    # pad number with leading zeros, for a total of 4 digits
    for i in range(len(numstr),4):
        numstr = '0' + numstr

    # construct a URL to be sent to the SDSS server (just doing eboss for now)
    url = 'https://dr14.sdss.org/sas/dr14/eboss/spectro/redux/v5_10_0/spectra/lite/'
    if online:
        url += plate + '/spec-' + plate + '-' + mjd + '-' + numstr + '.fits'

        # update the user what's going on
        if verbose:
            print("Trying to download SDSS spectrum from:\n {}\n".format(url))

        # send the URL to SDSS, and get back a response
        hdulist = pf.open(url)

    else:
        filein = 'spec-4055-55359-0596.fits'
        hdulist = pf.open(filein)

    # temporarily pull out wavelengths, fluxes, and the spectrum id
    logwav = hdulist[1].data['loglam']
    flux = hdulist[1].data['flux']
    ivar = hdulist[1].data['ivar']
    id = hdulist[2].data['SPECOBJID']

    # build spectrum dictionary
    s = {}

    # wavelength in angstroms
    s['wavelengths'] = np.array(10**logwav)

    # flux in 10e-17 erg cm^-2 s^-1 Angstrom^-1
    #  fluxconversion = 1.0e-3/1.0e-17/1e10
    s['flux'] = flux  # *fluxconversion

    # populate the metadata for this fake spectrum
    s['spectrumID'] = str(id[0])
    s['spectro_class'] = hdulist[2].data['class'][0]
    s['spectro_subclass'] = hdulist[2].data['subclass'][0]
    s['z'] = hdulist[2].data['z'][0]

    # which elements of the w and f arrays are in the optical?
    w = 10**logwav
    optical = (w > 3900) & (w < 8000)

    # is optical flux, on average brighter than 1.0e-17 erg/cm**2/s/angstrom?
    if np.mean(flux[optical]) > 1.0:
        # ...if so, return this spectrum
        if verbose:
            print("Returning a spectrum for {}.".format(id))
        return s
    else:
        # ...if not, call getSDSS() to try to get a brighter one
        if verbose:
            print("Yikes, {} was pretty faint. Trying again...".format(id))
        return getSDSS(online=online, verbose=verbose)


def getLabel(spectrum):
    '''
    Create a short string summarizing the basic properties of
    a SDSS spectrum. The spectrum must be in the dictionary
    format created by the SDSS web API.

        Arguments:

            spectrum = a dictionary containing an SDSS spectrum

        Returns:

            a string describing how the SDSS project has
            classified this object

        Example:
            # get a random spectrum
            s = getSDSS()

            # get its label, returned as a string
            label = getLabel(s)
            print label

    '''

    # create an empty label
    label = ""

    # add the identifier for this spectrum
    label += spectrum["spectrumID"] + "\n"

    # what is the basic class of this object?
    label += spectrum['spectro_class']

    # if there are more details, add them in
    if len(spectrum['spectro_subclass']) > 0:
        label += "/" + spectrum['spectro_subclass']
    label += ", "

    # show the redshift (z = v/c) of the object
    label += "redshift={:.2g}".format(spectrum["z"])

    # return that label as a string
    return label


h = 6.626068e-34    # Planck's constant, m**2 kg / s
k = 1.3806503e-23   # Boltzmann's constant, m**2 kg / s**2 / K
c = 299792458.0     # speed of light, m / s

def planck(w_nm, T=5800):
    """
    Compute the thermal emission spectrum.

    Arguments:
     w_nm = wavelength in nm, float or numpy array
     T = temperature in K (set to 5800K by default)

    Returns:
     intensity at the given wavelength(s), in W / m**3 / sr

    Example:
     >> import thermal
     >> print thermal.planck(550, T=5800)
     2.63106560099e+13
    """

    # convert nm to m
    w = w_nm/1.0e9

    # return Planck function
    num = 2 * h *c**2
    den = w**5 * (np.exp(h * c / (w * k * T)) - 1.0)
    return num/den

def fakeSDSS(T='random'):
    try:
        temperature = float(T)
    except (TypeError, ValueError):
        temperature = np.round(np.random.uniform(2000,9000))


    w = np.linspace(300, 1000, 1000)
    f = planck(w, T=temperature)

    # populate a fake SDSS spectrum
    s = {}

    # wavelength in angstroms
    s['wavelengths'] = np.array(w*10)

    # convert W/m^3 to
    fluxconversion = 1.0e-3/1.0e-17/1e10
    s['flux'] = f*fluxconversion

    # populate the metadata for this fake spectrum
    s['spectrumID'] = 'an imaginary SDSS'
    s['spectro_class'] = '{}K emission spectrum'.format(temperature)
    s['spectro_subclass'] = ''
    s['z'] = 0.0

    label = getLabel(s)
    return w, f, label

def testSDSS(**kwargs):
    '''This function tests the getSDSS() function
        and uses the fakeSDSS() function if that fails.'''

    # get a spectrum
    try:
        w, f, l = getSDSS(**kwargs)
        print("Hooray! It seems like the getSDSS() function works!")
    except:
        print("Uh-oh! It seems like something is broken between the SDSS server and you.")
        print(" Please use sdss.fakeSDSS() instead of sdss.getSDSS() to make spectra.")
        w, f, l = fakeSDSS()

    # plot a spectrum
    wmin = 390
    wmax = 800
    mask = (w>wmin) & (w<wmax)
    plt.plot(w[mask],f[mask])
    plt.title(l)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Flux ($10^{-17}$ erg cm$^{-2}$ s$^{-1}$ nm$^{-1}$)")
    plt.xlim(wmin, wmax)
    #plt.ylim(0, None)
    plt.show()

def getSDSS(**kwargs):
    '''
    Return a random spectrum from the Sloan Digital Sky Survey (SDSS).

    Returns
    -------
    wavelength : numpy array
        The wavelengths of the spectrum, in nm.
    flux : numpy array
        The fluxes of the spectrum, in 1e-17 erg/cm**2/s/nm
    description : string
        A string description explaining what the spectrum is.
    '''

    spectrum = downloadSDSS(**kwargs)
    label = getLabel(spectrum)

    wavelength_nm = spectrum['wavelengths']/10.0
    flux = spectrum['flux']*10

    return wavelength_nm, flux, label

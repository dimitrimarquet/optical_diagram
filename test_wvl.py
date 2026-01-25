import numpy as np
import matplotlib.pyplot as plt

def wavelength_to_rgb(wavelengths, gamma=0.8):
    '''This converts a given wavelength of light to an
    approximate RGB color value. The wavelength must be given
    in nanometers in the range from 380 nm through 750 nm
    (789 THz through 400 THz).
    Based on code by Dan Bruton
    http://www.physics.sfasu.edu/astro/color/spectra.html
    '''
    wavelengths = np.array(wavelengths).astype(float)
    rgbs = np.zeros(wavelengths.shape + (3,))

    attenuations = np.full(wavelengths.shape, np.nan)

    idx = ((wavelengths >= 380) & (wavelengths <= 440))
    attenuations[idx] = 0.3 + 0.7 * (wavelengths[idx] - 380) / (440 - 380)
    rgbs[idx,0] = ((-(wavelengths[idx] - 440) / (440 - 380)) * attenuations[idx]) ** gamma
    rgbs[idx,2] = (1.0 * attenuations[idx]) ** gamma

    idx = ((wavelengths >= 440) & (wavelengths <= 490))
    rgbs[idx,1] = ((wavelengths[idx] - 440) / (490 - 440)) ** gamma
    rgbs[idx,2] = 1.0

    idx = ((wavelengths >= 490) & (wavelengths <= 510))
    rgbs[idx,1] = 1.0
    rgbs[idx,2] = (-(wavelengths[idx] - 510) / (510 - 490)) ** gamma

    idx = ((wavelengths >= 510) & (wavelengths <= 580))
    rgbs[idx,0] = ((wavelengths[idx] - 510) / (580 - 510)) ** gamma
    rgbs[idx,1] = 1.0

    idx = ((wavelengths >= 580) & (wavelengths <= 645))
    rgbs[idx,0] = 1.0
    rgbs[idx,1] = (-(wavelengths[idx] - 645) / (645 - 580)) ** gamma

    idx = ((wavelengths >= 645) & (wavelengths <= 750))
    attenuations[idx] = 0.3 + 0.7 * (750 - wavelengths[idx]) / (750 - 645)
    rgbs[idx,0] = (1.0 * attenuations[idx]) ** gamma

    rgbs = (255*rgbs).astype(int)
    return rgbs

print(wavelength_to_rgb(700))

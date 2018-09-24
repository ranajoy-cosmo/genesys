#!/usr/bin/env python

import numpy as np
import sys

def get_maj_min_fwhm(fwhm_mid, elp):
    x = ((2 + elp)/(2 - elp))**2
    maj = np.sqrt(2)*fwhm_mid*x/np.sqrt(1 + x)
    min = np.sqrt(2)*fwhm_mid/np.sqrt(1 + x)
    return maj, min

if __name__=="__main__":
    mid = np.float(sys.argv[1])
    elp = np.float(sys.argv[2])
    fwhm_maj, fwhm_min = get_maj_min_fwhm(mid, elp)
    print fwhm_maj, fwhm_min

#!/usr/bin/env python

import numpy as np
import healpy as hp

pointing_p_params = {'alpha' : 45.0, 'beta' : 45.0, 'sampling_rate' : 10.0, 'segment_length' : 3600.0, 't_spin' : 60.0, 't_precession' : 4*24*60*60.0, 't_year' : 365.25*24*60*60.0, 'fp_pos' : [0.0,0.0], 'offset' : [0.0,0.0], 'coord_system' : 'ecliptic', 'pol_phase_ini' : 0.0}

config = Generic_Class()
config.__dict__.update(config_dict)

pt_obj = pt.Pointing(config, 0, 1)
vec_obv = pt_obj.get_vec_obv(0.0)
pol_ang = pt_obj.get_pol_ang(vec_obv)

def simple_binning(vec_obv, nside=128):
    npix = 12*nside**2
    hitpix = hp.vec2pix(nside, vec_obv[...,0], vec_obv[...,1], vec_obv[...,2])
    hitmap = np.bincount(hitpix, minlength=npix)
    return hitmap

hitmap_2 = simple_binning(vec_obv, 128)

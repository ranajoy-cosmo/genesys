#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
import importlib
from simulation.lib.plotting.my_imshow import new_imshow
from scipy.signal import convolve2d
import convolution_kernel

def gaussian_2d(beam_params, bolo_params, mesh):
    factor = 2*np.sqrt(2*np.log(2))                                     #FWHM <--> sigma
    sigma = bolo_params.fwhm/factor

    #Building a circular Gaussian beam on the 2D mesh
    x,y = mesh
    beam_kernel = np.exp(-(x**2/(2*sigma**2) + y**2/(2*sigma**2)))
    
    #Normalising the circular Gaussian beam so that the integral is 1
    integral = np.sum(beam_kernel)*beam_params.beam_resolution**2
    #integral = 2*np.pi*sigma**2
    beam_kernel /= integral

    if bolo_params.conv_fwhm == 0.0:
        return beam_kernel, None
        #return <beam_kernel>, <conv_kernel>
    else:
        conv_kernel, del_x = convolution_kernel.get_beam(beam_params, bolo_params)
        beam_kernel = convolve2d(beam_kernel, conv_kernel, mode="same")*beam_params.beam_resolution
    return beam_kernel, conv_kernel


def check_normalisation(beam_params, bolo_params, beam_kernel):
    dx = beam_params.beam_resolution
    dy = beam_params.beam_resolution
    integral = np.sum(beam_kernel)*dx*dy
    print "The integral of the beam is :", integral
    print "Percentage difference with unity :", 100*(1-integral)


def get_mesh(beam_params, bolo_params):
    offset_max = max(abs(bolo_params.del_x), abs(bolo_params.del_y))              #arc-mins
    fwhm_minor = bolo_params.fwhm
    fwhm_major = np.sqrt(fwhm_minor**2 + bolo_params.conv_fwhm**2)
    size = beam_params.beam_cutoff*fwhm_major + offset_max             #arc-mins
    dd = beam_params.beam_resolution                                            #arc-mins
    n = int(size/dd/2)
    x = np.arange(-n, n+1)*dd
    y = -1*np.arange(-n, n+1)*dd
    return np.meshgrid(x,y), x 


def display_beam_settings(beam_params, bolo_params, mesh):
    if beam_params.do_pencil_beam:
        print "Pencil beam"
    else:
        factor = 2*np.sqrt(2*np.log(2))
        fwhm_minor = bolo_params.fwhm
        fwhm_major = np.sqrt(fwhm_minor**2 + bolo_params.conv_fwhm**2)
        ellipticity = 100*2*(fwhm_major - fwhm_minor)/(fwhm_major + fwhm_minor)
        print "Major axis(FWHM) :", fwhm_major, "arcmins" 
        print "Minor axis(FWHM) :", fwhm_minor, "arcmins"
        print "Ellipticity :", ellipticity, "%"
        print "Center :", bolo_params.del_x, bolo_params.del_y
        print "Tilt :", bolo_params.beam_angle, "degrees"
        print "Pixel size :", beam_params.beam_resolution, "arcmins" 
        print "Kernel width in FWHM of beam:", beam_params.beam_cutoff
        print "# of pixels per FWHM (minor-axis) of beam :", fwhm_minor/beam_params.beam_resolution
        print "Expected # of pixels in kernel cross-section :", int(beam_params.beam_cutoff*fwhm_major/beam_params.beam_resolution/2)*2 + 1 
        print "Actual # of pixels in kernel cross-section :", mesh[0][0].size 

def plot_beam(beam_kernel, beam_params):
    fig, ax = plt.subplots()
    n = beam_kernel[0].size/2
    extent = np.arange(-n, n+1)*beam_params.beam_resolution
    im = new_imshow(ax, beam_kernel, x=extent, y=extent, interpolation="nearest")
    fig.colorbar(im, ax=ax)
    plt.show()


if __name__=="__main__":

    from custom_params import beam_params
    bolo_name = sys.argv[1]
    bolo_params = importlib.import_module("simulation.timestream_simulation.bolo_params." + bolo_name).bolo_params

    if beam_params.do_pencil_beam:
        beam_kernel = np.array([[1]])/beam_params.beam_resolution**2
        del_beta = np.array([0])
        beam_params.plot_beam = False
    else:
        mesh, del_beta = get_mesh(beam_params, bolo_params)
        beam_kernel, convolve_kernel = gaussian_2d(beam_params, bolo_params, mesh)

    if beam_params.check_normalisation:
        check_normalisation(beam_params, bolo_params, beam_kernel)
    if beam_params.display_beam_settings:
        display_beam_settings(beam_params, bolo_params, mesh)
    if beam_params.plot_beam:
        plot_beam(beam_kernel, beam_params)
        if bolo_params.ellipticity != 0.0:
            plot_beam(convolve_kernel, beam_params)
    print del_beta


def get_beam(beam_params, bolo_params):
    if beam_params.do_pencil_beam:
        beam_kernel = np.array([[1]])
        del_beta = np.array([0])
    else:
        mesh, del_beta = get_mesh(beam_params, bolo_params)
        beam_kernel, convolve_kernel = gaussian_2d(beam_params, bolo_params, mesh)
    return beam_kernel, del_beta

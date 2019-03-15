#!/usr/bin/env python

import healpy as hp
import matplotlib.pyplot as plt
import os
import sys
from genesys.plotting.colourmaps import cmap_planck_cmb

def save_map(sky_map, output_dir, output_file_name, range_min=None, range_max=None, unit="$\mu$K", cmap=cmap_planck_cmb, pol=True, sky_mask=None):

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    figsize_inch = 12, 8
    dpi = 100

    if sky_mask is not None:
        sky_map *= sky_mask

    if pol == True:
        if range_min == None:
            range_min = (None, None, None)
        if range_max == None:
            range_max = (None, None, None)
        for (label, num, r_min, r_max) in zip('IQU', (0,1,2), range_min, range_max):
            fig = plt.figure(figsize=figsize_inch, dpi=dpi)
            hp.mollview(sky_map[num], fig=fig.number, min=r_min, max=r_max, unit=unit, cmap=cmap, xsize=figsize_inch[0]*dpi, title=output_file_name+'_'+label)
            plt.savefig(os.path.join(output_dir, output_file_name+'_'+label+'.png'), dpi=dpi, bbox_inches="tight")
    else:
        fig = plt.figure(figsize=figsize_inch, dpi=dpi)
        hp.mollview(sky_map, fig=fig.number, min=range_min, max=range_max, unit=unit, cmap=cmap, xsize=figsize_inch[0]*dpi, title=output_file_name)
        plt.savefig(os.path.join(output_dir, output_file_name+'.png'), dpi=dpi, bbox_inches="tight")

if __name__=="__main__":
    sky_map_name = sys.argv[1]
    mask_name = sys.argv[2]

    try:
        sky_map = hp.read_map(sky_map_name, field=(0,1,2), verbose=False)
    except IndexError:
        sky_map = hp.read_map(sky_map_name, verbose=False)

    try:
        mask = hp.read_map(mask_name, field=(0,1,2), verbose=False)
    except IndexError:
        mask = hp.read_map(mask_name, verbose=False)

    plot_map(sky_map, sky_mask=sky_mask, range_min=None, range_max=None, pol=False)

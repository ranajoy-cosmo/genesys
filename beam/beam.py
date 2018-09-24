#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.ndimage import zoom, interpolation
import sys
import os
import importlib
from genesys.utilities.plotting.my_imshow import new_imshow
from genesys.utilities import Generic_Class

class Beam():
    def __init__(self, config, bolo_config):
        self.config = Generic()
        self.config.__dict__.update(config.__dict__)
        self.config.__dict__.update(bolo_config.__dict__)
        self.get_beam()


    def get_beam(self):
        if self.config.beam_type == "pencil":
            self.beam_kernel = np.array([[1.0]])
            self.del_beta = np.array([0.0])
        elif self.config.beam_type == "full_simulated":
            mesh = self.get_mesh()
            self.gaussian_2d(mesh)
        elif self.config.beam_type == "from_file":
            self.read_mark_beam_map()
        else:
            print("Beam type not recognised.")
            sys.exit()


    def read_mark_beam_map(self):
        beam_kernel = np.load(self.config.input_beam_file)
        mark_orig_dim = 201                             #pixels
        mark_fwhm_major = 7.58                   #arc-mins
        mark_fwhm_minor = 7.58                   #arc-mins
        mark_resolution = 0.30                   #arc-mins

        new_dim = int(mark_orig_dim * mark_resolution / self.config.scan_resolution)
        if new_dim%2 == 0:
            new_dim += 1
        #print "old resolution :", mark_resolution
        #print "new_resolution :", self.config.scan_resolution
        old_extent = mark_orig_dim * mark_resolution
        new_extent = new_dim * self.config.scan_resolution
        #print "old extent :", mark_resolution * 181
        #print "new_extent :", new_extent
        self.config.fwhm_major = mark_fwhm_major * new_extent / old_extent
        self.config.fwhm_minor = mark_fwhm_minor * new_extent / old_extent
        #print "old fwhm :", mark_fwhm_major
        #print "new fwhm :", self.config.fwhm_major

        self.beam_kernel = np.empty((4, new_dim, new_dim))
        for i in range(4):
            self.beam_kernel[i] = zoom(beam_kernel[i], float(new_dim)/float(mark_orig_dim))
            self.beam_kernel[i] = interpolation.rotate(self.beam_kernel[i], angle=self.config.beam_angle, reshape=False)

        #print "beam_cutoff required :", self.config.beam_cutoff
        beam_extension_required = self.config.fwhm_major * self.config.beam_cutoff          #arc-min
        #print "extention required :", beam_extension_required
        num_pix = int(new_dim * beam_extension_required / new_extent)
        if num_pix%2 == 0:
            num_pix -= 1
        #print "num pix :", num_pix
        #print "extension achieved :", num_pix * self.config.scan_resolution

        self.config.beam_cutoff = num_pix * self.config.scan_resolution / self.config.fwhm_major
        #print "cutoff achieved :", self.config.beam_cutoff

        start = new_dim/2 - num_pix/2
        stop = new_dim/2 + num_pix/2 + 1
        self.beam_kernel = self.beam_kernel[:, start:stop, start:stop]
        #print "kernel shape :", self.beam_kernel.shape
        self.del_beta = self.config.scan_resolution * np.arange(-num_pix/2 + 1, num_pix/2 + 1)


    def gaussian_2d(self, mesh):
        factor = 2*np.sqrt(2*np.log(2))                                     #FWHM <--> sigma

        sigma_major = self.config.fwhm_major/factor
        sigma_minor = self.config.fwhm_minor/factor
        #theta = np.deg2rad(self.config.beam_angle + self.config.pol_phase_ini)
        theta = np.deg2rad(self.config.beam_angle)
        x0, y0 = self.config.offset_x/60.0, self.config.offset_y/60.0

        #Building a circular Gaussian beam on the 2D mesh
        x,y = mesh
        a = (np.cos(theta)**2)/(2*sigma_major**2) + (np.sin(theta)**2)/(2*sigma_minor**2)
        b = 1*np.sin(2*theta)/(4*sigma_major**2) - np.sin(2*theta)/(4*sigma_minor**2)
        c = (np.sin(theta)**2)/(2*sigma_major**2) + (np.cos(theta)**2)/(2*sigma_minor**2)
        norm_factor = 2*np.pi*sigma_major*sigma_minor
        beam_kernel = np.exp(-1*(a*(x - x0)**2 + 2*b*(x - x0)*(y - y0) + c*(y - y0)**2)) 
        beam_kernel /= norm_factor
        
        size = beam_kernel.shape[0]
        self.beam_kernel = np.zeros((4, size, size))
        self.beam_kernel[0] = beam_kernel
        self.beam_kernel[1] = beam_kernel


    def get_beam_row(self, del_beta):
        row_num = np.where(self.del_beta==del_beta)[0][0]

        return self.beam_kernel[:,row_num]

    def block_beam_row(self, del_beta):
        row_num = np.where(self.del_beta==del_beta)[0][0]

        self.beam_kernel[:,row_num] = np.zeros(self.del_beta.size)


    def normalise(self):
        dx = self.config.scan_resolution
        dy = self.config.scan_resolution
        integral = np.sum(self.beam_kernel[0])*dx*dy
        for i in range(4): 
            self.beam_kernel[i] /= integral


    def check_normalisation(self):
        dx = self.config.scan_resolution
        dy = self.config.scan_resolution
        map_index = ["T", "Q", "U", "V"]
        for i in range(4): 
            integral = np.sum(self.beam_kernel[i])*dx*dy
            print("The integral of the {} beam is : {}".format(map_index[i], integral))
            print("Percentage difference with unity : {}".format(100*(1-integral)))


    def get_mesh(self):
        offset_max = max(abs(self.config.offset_x/60.0), abs(self.config.offset_y/60.0))              #arc-mins
        fwhm_major = self.config.fwhm_major
        fwhm_minor = self.config.fwhm_minor
        size = self.config.beam_cutoff*fwhm_major + offset_max             #arc-mins
        dd = self.config.scan_resolution                                            #arc-mins
        n = int(size/dd/2)
        x = np.arange(-n, n+1)*dd
        y = -1*np.arange(-n, n+1)*dd
        self.del_beta = y
        return np.meshgrid(x,y)

    def shift_beam_kernel(self, direction, units=1, fill_value=0):
        if direction == "top":
            self.old_beam_kernel = np.copy(self.beam_kernel) 
            self.beam_kernel = np.full(self.beam_kernel.shape, fill_value)
            beam_dim = self.del_beta.size
            for i in range(4):
                self.beam_kernel[i, :beam_dim-units] = self.old_beam_kernel[i, units:] 

        if direction == "bottom":
            self.old_beam_kernel = np.copy(self.beam_kernel) 
            self.beam_kernel = np.full(self.beam_kernel.shape, fill_value)
            beam_dim = self.del_beta.size
            for i in range(4):
                self.beam_kernel[i, units:] = self.old_beam_kernel[i, :beam_dim-units] 

    def restore_beam_kernel(self):
        self.beam_kernel = np.copy(self.old_beam_kernel)

    def display_beam_settings(self):
        factor = 2*np.sqrt(2*np.log(2))
        fwhm_major = self.config.fwhm_major
        fwhm_minor = self.config.fwhm_minor
        ellipticity = 100*2*(fwhm_major - fwhm_minor)/(fwhm_major + fwhm_minor)
        display_string = "\n#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
        display_string += "#* BEAM PARAMETERS : {}\n".format(self.config.name) 
        display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n"
        display_string += "Major axis(FWHM) : {} arcmins\n".format(fwhm_major) 
        display_string += "Minor axis(FWHM) : {} arcmins\n".format(fwhm_minor)
        display_string += "Ellipticity : {} %\n".format(ellipticity)
        display_string += "Center : ({},{})\n".format(self.config.offset_x, self.config.offset_y)
        display_string += "Tilt : {} degrees\n".format(self.config.beam_angle)
        display_string += "Kernel width in FWHM of beam : {}\n".format(self.config.beam_cutoff)
        display_string += "# of pixels per FWHM (minor-axis) of beam : {}\n".format(fwhm_minor/self.config.scan_resolution)
        display_string += "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*\n\n"
        prompt(display_string, sys.stdout)

    """
    def plot_beam(self):
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        n = self.beam_kernel[0].shape[0]/2
        extent = np.arange(-n, n+1)*self.config.scan_resolution
        im = new_imshow(ax1, self.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
        ax1.set_title('T')
        fig.colorbar(im, ax=ax1)
        im = new_imshow(ax2, self.beam_kernel[1], x=extent, y=extent, interpolation="nearest")
        ax2.set_title('Q')
        fig.colorbar(im, ax=ax2)
        im = new_imshow(ax3, self.beam_kernel[2], x=extent, y=extent, interpolation="nearest")
        ax3.set_title('U')
        fig.colorbar(im, ax=ax3)
        im = new_imshow(ax4, self.beam_kernel[3], x=extent, y=extent, interpolation="nearest")
        ax4.set_title('V')
        fig.colorbar(im, ax=ax4)
        fig.suptitle("Rescaled Plack 217_5a, FWHM : 7.68', Resol : 0.96', Extent : 3.85*FWHM")
        plt.show()


    def plot_beam(self):
        fig, ax= plt.subplots()
        n = self.beam_kernel[0].shape[0]/2
        extent = np.arange(-n, n+1)*self.config.scan_resolution
        im = new_imshow(ax, 10*self.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
        fig.colorbar(im)
        #im = new_imshow(ax, 10*self.beam_kernel[1], x=extent, y=extent, interpolation="nearest", cmap='gray')
        plt.show()
    """

    def write_beam(self, out_dir=None):
        if out_dir == None:
            out_dir = os.path.join(os.getcwd(), "beam_maps") 
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        np.save(os.path.join(out_dir, self.config.beam_file_name), self.beam_kernel)


def plot_beam(bolo_1, bolo_2, bolo_3):
    fig, ((ax11, ax12, ax13), (ax21, ax22, ax23), (ax31, ax32, ax33)) = plt.subplots(3, 3, sharex='col', sharey='row')
    n = bolo_1.beam_kernel[0].shape[0]/2
    extent = np.arange(-n, n+1)*bolo_1.config.scan_resolution

    im = new_imshow(ax11, bolo_1.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
    ax11.set_title('Central')
    #fig.colorbar(im, ax=ax11)
    im = new_imshow(ax12, bolo_2.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
    ax12.set_title('Top')
    #fig.colorbar(im, ax=ax12)
    im = new_imshow(ax13, bolo_3.beam_kernel[0], x=extent, y=extent, interpolation="nearest")
    ax13.set_title('Bottom')
    fig.colorbar(im, ax=ax13)

    im = new_imshow(ax21, bolo_1.beam_kernel[1], x=extent, y=extent, interpolation="nearest")
    #ax21.set_title('T')
    #fig.colorbar(im, ax=ax21)
    im = new_imshow(ax22, bolo_2.beam_kernel[1], x=extent, y=extent, interpolation="nearest")
    #ax22.set_title('Q')
    #fig.colorbar(im, ax=ax22)
    im = new_imshow(ax23, bolo_3.beam_kernel[1], x=extent, y=extent, interpolation="nearest")
    #ax23.set_title('U')
    fig.colorbar(im, ax=ax23)

    im = new_imshow(ax31, 1000*bolo_1.beam_kernel[2], x=extent, y=extent, interpolation="nearest")
    #ax31.set_title('T')
    #fig.colorbar(im, ax=ax31)
    im = new_imshow(ax32, 1000*bolo_2.beam_kernel[2], x=extent, y=extent, interpolation="nearest")
    #ax32.set_title('Q')
    #fig.colorbar(im, ax=ax32)
    im = new_imshow(ax33, 1000*bolo_3.beam_kernel[2], x=extent, y=extent, interpolation="nearest")
    #ax33.set_title('U')
    fig.colorbar(im, ax=ax33)

    fig.text(0.04, 0.78, 'T', ha='center', va='center')
    fig.text(0.04, 0.5, 'Q', ha='center', va='center')
    fig.text(0.05, 0.22, r'U$\times$1000', ha='center', va='center')
    fig.text(0.075, 0.5, 'arcmins', ha='center', va='center', rotation='vertical')
    fig.text(0.5, 0.05, 'arcmins', ha='center', va='center')

    fig.suptitle("Simulated CORE beams using GRASP")
    plt.rc('xtick', labelsize=8)
    plt.rc('ytick', labelsize=8)
    plt.show()

if __name__=="__main__":

    config_file = sys.argv[1]
    #bolo_name = sys.argv[2]
    config = importlib.import_module("simulation.timestream_simulation.config_files." + config_file).config
    config.scan_resolution = 0.3 
    #bolo_config = importlib.import_module("simulation.timestream_simulation.bolo_config_files." + config.bolo_config_file).bolo_config
    from simulation.timestream_simulation.bolo_config_files.four_bolos_optimal import bolo_config

    bolo_beam_1 = Beam(config, bolo_config.bolos[config.bolo_list[0]])
    bolo_beam_2 = Beam(config, bolo_config.bolos[config.bolo_list[1]])
    bolo_beam_3 = Beam(config, bolo_config.bolos[config.bolo_list[2]])

    #if config.check_normalisation:
    #    bolo_beam.check_normalisation()
    #if config.display_beam_settings:
    #    bolo_beam.display_beam_settings()
    #bolo_beam.plot_beam()
    bolo_beam_1.normalise()
    #bolo_beam_1.check_normalisation()
    bolo_beam_2.normalise()
    #bolo_beam_2.check_normalisation()
    bolo_beam_3.normalise()
    #bolo_beam_3.check_normalisation()

    plot_beam(bolo_beam_1, bolo_beam_2, bolo_beam_3)
    #if config.write_beam:
    #    bolo_beam.write_beam()

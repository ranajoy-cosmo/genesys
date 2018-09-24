import numpy as np
import healpy as hp
import sys
from genesys.utilities import Generic_Class
from genesys.utilities import prompt
import genesys.timestream_simulation.pointing.quaternion.quaternion as qt

class Pointing():
    
    def __init__(self, config, segment, beam_width):
        self.config = Generic()
        self.config.__dict__.update(config.__dict__)
        self.beam_width = beam_width
        self.segment = segment
        self.set_n_samples()
        self.set_sim_times()
        self.set_initial_axes_and_angles()
        self.generate_quaternion()

    def set_n_samples(self):
        if self.config.beam_type == "pencil_beam":
            self.pad = 0
        else:
            self.pad = self.beam_width/2 

        if self.config.sim_type == "template":
            if self.config.template_type == "tm_gradient":
                self.pad += 1

#        prompt("{}, {}, {}, {}\n".format(int(self.config.t_segment*self.config.sampling_rate), int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate, int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate + 2*self.pad, self.pad))
        self.nsamples = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate + 2*self.pad

    def set_sim_times(self):
        self.t_start = self.segment*self.config.t_segment
        self.t_stop = (self.segment + 1)*self.config.t_segment
        self.delta_t = 1.0/self.config.sampling_rate/self.config.oversampling_rate

    #Setting the initial axes
    def set_initial_axes_and_angles(self):
        self.alpha = np.radians(self.config.alpha)                                   #radians
        self.beta = np.radians(self.config.beta)                                     #radians

        self.centre_opening = self.alpha + self.beta + self.config.focal_plane_del_beta

        self.axis_spin = np.array([np.cos(self.alpha), 0.0, np.sin(self.alpha)])
        self.axis_prec = np.array([1.0, 0.0, 0.0])
        self.axis_rev = np.array([0.0, 0.0, 1.0])


    def generate_quaternion(self):
        t_steps = self.t_start + self.delta_t*np.arange(-self.pad, self.nsamples - self.pad)

        w_spin = 2*np.pi/self.config.t_spin
        w_prec = 2*np.pi/self.config.t_prec
        w_rev = 2*np.pi/self.config.t_year

        self.r_total = quaternion.multiply(quaternion.make_quaternion(w_rev*t_steps, self.axis_rev), quaternion.multiply(quaternion.make_quaternion(w_prec*t_steps, self.axis_prec), quaternion.make_quaternion(w_spin*t_steps, self.axis_spin)))


    #Setting the initial positioning of pointing vectors
    def get_initial_vec(self, del_beta):
        del_beta_rad = np.deg2rad(del_beta/60.0)                                #radians

        if self.config.beam_type == "pencil":
            del_x = np.radians(self.config.offset_x/60.0/60.0)    #radians
            del_y = np.radians(self.config.offset_y/60.0/60.0)    #radians
        else:
            del_x = 0.0
            del_y = 0.0

        total_opening = self.centre_opening + del_beta_rad

        u_view = np.array([np.cos(total_opening), 0.0, np.sin(total_opening)])
        
        x_roll_axis = np.array([0.0, 1.0, 0.0])
        y_roll_axis = np.array([-np.sin(total_opening), 0.0, np.cos(total_opening)])

        q_x_roll = quaternion.make_quaternion(del_x, x_roll_axis)
        q_y_roll = quaternion.make_quaternion(del_y, y_roll_axis)
        q_offset = quaternion.multiply(q_x_roll, q_y_roll)

        u_view = quaternion.transform(q_offset, u_view)

        return u_view

    #Simulating the pointing
    def get_vec_obv(self, del_beta):
        v_init = self.get_initial_vec(del_beta)
        v = quaternion.transform(self.r_total, v_init)

        if self.config.coordinate_system == "galactic":
            v = self.transform_to_gal_coords(v)

        return v

    def transform_to_gal_coords(self, v):
        rot = hp.Rotator(coord=['E', 'G'])
        theta, phi = hp.vec2ang(v)
        theta_gal, phi_gal = rot(theta, phi)
        del theta, phi
        v = hp.ang2vec(theta_gal, phi_gal)
        del theta_gal, phi_gal
        return v 


    def get_pol_ang(self, v_pointing):
        pol_ini = np.radians(self.config.pol_phase_ini)
        pol_vec_ini = np.array([0.0, 1.0, 0.0])

        pol_vec = quaternion.transform(self.r_total, np.tile(pol_vec_ini, self.nsamples).reshape(-1,3))
        if self.config.coordinate_system == "galactic":
            pol_vec = self.transform_to_gal_coords(pol_vec)

        theta, phi = hp.vec2ang(v_pointing)

        x_local = np.array(zip(np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)))
        y_local = np.array(zip(-np.sin(phi), np.cos(phi), np.zeros(phi.size)))

        proj_x = np.sum(pol_vec*x_local, axis=-1)
        proj_y = np.sum(pol_vec*y_local, axis=-1)
        del pol_vec, x_local, y_local

        #pol_ang = np.pi - (np.arctan2(proj_y, proj_x) + pol_ini) % np.pi 
        pol_ang = (np.arctan2(proj_y, proj_x) + pol_ini) % np.pi 
        del proj_x, proj_y

        return pol_ang 

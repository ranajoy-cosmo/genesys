import numpy as np
import healpy as hp
import sys
from . import pointing_utils as pu
from ...utilities import Generic_Class, prompt

class Pointing():
    
    def __init__(self, config, segment, beam_width):
        self.config = Generic_Class()
        self.config.__dict__.update(config.__dict__)
        self.beam_width = beam_width
        self.segment = segment
        self.set_n_samples()
        self.set_sim_times()
        self.set_initial_axes_and_angles()
        self.generate_quaternion()

    def set_n_samples(self):
        """
        Set the number of samples in the timestream.
        If an extended pixelised beam is used, a pad of length (beam_width - 1)/2 is added on either side.
        The beam width will always be an odd integer.
        If a gradient is to be done, an additional element is added for the padding.
        """
        if self.config.beam_type == "pencil_beam":
            self.pad = 0
        else:
            self.pad = int(self.beam_width / 2)

        if self.config.sim_type == "template" and self.config.template_type == "tm_gradient":
            self.pad += 1

        # The time quantities are in seconds and the sampliong quantities are in Hz
        self.nsamples = int(self.config.t_segment*self.config.sampling_rate)*self.config.oversampling_rate + 2*self.pad

    def set_sim_times(self):
        """
        Set the start and end times for the given segment.
        Also set the time interval between samples.
        All quantities are in seconds
        """
        self.t_start = self.segment*self.config.t_segment
        self.t_stop = (self.segment + 1)*self.config.t_segment
        self.delta_t = 1.0 / self.config.sampling_rate / self.config.oversampling_rate

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    def set_initial_axes_and_angles(self):
        """
        The initial position of the satellite is assumed to be such that the
        anti-solar precession axis is the x-axis,
        the spin axis is in the x-z plane and at a positive angle of alpha degrees to the precession axis,
        the axis of revolution is the z-axis.
        The initial pointing of individual detectors is set in gen_initial_vec()
        """
        # alpha and beta are provided in degrees in the config files.
        # They have to be converted to radians
        alpha = np.radians(self.config.alpha)                                   #radians
        beta = np.radians(self.config.beta)                                     #radians
        boresight_opening_angle = alpha + beta

        self.axis_rev = np.array([0.0, 0.0, 1.0])
        self.axis_prec = np.array([1.0, 0.0, 0.0])
        self.axis_spin = np.array([np.cos(alpha), 0.0, np.sin(alpha)])
        self.axis_boresight = np.array([np.cos(boresight_opening_angle), 0.0, np.sin(boresight_opening_angle)])
        self.axis_boresight_rot_x = np.array([np.cos(boresight_opening_angle - np.pi/2.0), 0.0, np.sin(boresight_opening_angle - np.pi/2.0)])

    def generate_quaternion(self):
        t_steps = self.t_start + self.delta_t*np.arange(-self.pad, self.nsamples - self.pad)

        w_rev = 2*np.pi/self.config.t_year
        w_prec = 2*np.pi/self.config.t_prec
        w_spin = 2*np.pi/self.config.t_spin

        self.total_rotation_quaternion = pu.gen_rotation_quat(w_rev*t_steps, self.axis_rev) * pu.gen_rotation_quat(w_prec*t_steps, self.axis_prec) * pu.gen_rotation_quat(w_spin*t_steps, self.axis_spin)

        if self.config.coord_system == "galactic":
            self.total_rotation_quaternion = pu.quaternion_coordinate_transformation('E','G') * self.total_rotation_quaternion

    def get_initial_vec(self, beam_row_del_beta=0):
        """
        The initial pointing of each detector is set independently.
        The position on the focal plane for each detector is different, hence their projection on the sky.
        """
        # beam_row_del_beta is non-zero only when using an extended, real-space beam.
        # The direction of the beam_row_del_beta is same as the direction of beta.
        beam_row_del_beta_rad = np.radians(beam_row_del_beta / 60.0)                                #radians

        # The focal plane position, fp_pos is provided in the config in arc-mins
        fp_pos_x_rad = np.radians(self.config.fp_pos[0] / 60.0)
        fp_pos_y_rad = np.radians(self.config.fp_pos[1] / 60.0)
        # The offsets are given in arc-secs
        offset_x_rad = np.radians(self.config.offset[0] / 60.0 / 60.0)
        offset_y_rad = np.radians(self.config.offset[1] / 60.0 / 60.0)

        total_fp_pos_x_rad = fp_pos_x_rad + offset_x_rad
        total_fp_pos_y_rad = fp_pos_y_rad + offset_y_rad + beam_row_del_beta_rad

        y_axis = np.array([0.0,1.0,0.0])
        q_rot_x = pu.gen_rotation_quat(total_fp_pos_x_rad, self.axis_boresight_rot_x)
        q_rot_y = pu.gen_rotation_quat(total_fp_pos_y_rad, -1.0*y_axis)

        q_rot_total = q_rot_x * q_rot_y

        pointing_vector_initial = pu.rotate_vector(q_rot_total, self.axis_boresight)

        return pointing_vector_initial

    #Simulating the pointing
    def get_vec_obv(self, beam_row_del_beta):
        v_init = self.get_initial_vec(beam_row_del_beta)
        v_observed = pu.rotate_vector(self.total_rotation_quaternion, v_init)
        return v_observed


    def get_pol_ang(self, v_pointing):
        pol_ini = np.radians(self.config.pol_phase_ini)
        pol_vec_ini = np.array([0.0, 1.0, 0.0])

        pol_vec = pu.rotate_vector(self.total_rotation_quaternion, pol_vec_ini)

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

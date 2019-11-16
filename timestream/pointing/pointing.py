import sys
import math
import numpy as np
import healpy as hp
from . import pointing_utils as pu
from genesys import Genesys_Class

class Pointing(Genesys_Class):
    
    def __init__(self, params):
        self.params.__dict__.update(params.__dict__)
        self.set_initial_axes_and_angles()
        self.set_initial_pointing_vector()

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

    def set_initial_pointing_vector(self, beam_row_del_beta=0):
        """
        The initial pointing of each detector is set independently.
        The position on the focal plane for each detector is different, hence their projection on the sky.
        """
        # beam_row_del_beta is non-zero only when using an extended, real-space beam.
        # The direction of the beam_row_del_beta is same as the direction of beta.
        beam_row_del_beta_rad = np.radians(beam_row_del_beta / 60.0)                                #radians

        # The focal plane position, fp_pos is provided in the config in arc-mins
        fp_pos_x_rad = np.radians(self.config.focal_plane_pos[0] / 60.0)
        fp_pos_y_rad = np.radians(self.config.focal_plane_pos[1] / 60.0)
        # The offsets are given in arc-secs
        offset_x_rad = np.radians(self.config.offset[0] / 60.0 / 60.0)
        offset_y_rad = np.radians(self.config.offset[1] / 60.0 / 60.0)

        total_fp_pos_x_rad = fp_pos_x_rad + offset_x_rad
        total_fp_pos_y_rad = fp_pos_y_rad + offset_y_rad + beam_row_del_beta_rad

        y_axis = np.array([0.0,1.0,0.0])
        q_rot_x = pu.gen_rotation_quat(total_fp_pos_x_rad, self.axis_boresight_rot_x)
        q_rot_y = pu.gen_rotation_quat(total_fp_pos_y_rad, -1.0*y_axis)

        q_rot_total = q_rot_x * q_rot_y

        self.pointing_vector_initial = pu.rotate_vector(q_rot_total, self.axis_boresight)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def initialise_for_segment(self, segment):
        """
        This is called for each new segment
        The simulation times, quaternions are set accordingly for the particular segment
        """
        self.segment = segment
        self.generate_quaternion()

    def get_t_steps(self):
        delta_t = 1.0 / self.config.sampling_rate
        t_start = self.segment[1]
        t_stop = self.segment[2]
        t_duration = t_stop - t_start
        self.n_samples = int(math.ceil(t_duration / delta_t))
        t_steps = t_start + delta_t * np.arange(self.n_samples)
        return t_steps

    def generate_quaternion(self):
        t_steps = self.get_t_steps()

        w_rev = 2*np.pi/self.config.t_year
        w_prec = 2*np.pi/self.config.t_prec
        w_spin = 2*np.pi/self.config.t_spin

        self.total_rotation_quaternion = pu.gen_rotation_quat(w_rev*t_steps, self.axis_rev) * pu.gen_rotation_quat(w_prec*t_steps, self.axis_prec) * pu.gen_rotation_quat(w_spin*t_steps, self.axis_spin)

        if self.config.coordinate_system == "galactic":
            self.total_rotation_quaternion = pu.quaternion_coordinate_transformation('E','G') * self.total_rotation_quaternion

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    # Simulating the pointing
    def get_vec_obv(self):
        v_observed = pu.rotate_vector(self.total_rotation_quaternion, self.pointing_vector_initial)
        return v_observed

    # Generating the polarisation angle
    def get_pol_ang(self, v_observed):
        if self.config.polarisation_modulation == "instrumental_scanning":
            initial_phase = np.radians(self.config.pol_phase_ini)
            y_axis = pu.rotate_vector(self.total_rotation_quaternion, np.array([0.0,1.0,0.0]))
            sinalpha = y_axis[...,0] * v_observed[...,1] - y_axis[...,1] * v_observed[...,0]
            cosalpha = y_axis[...,2] - v_observed[...,2] * np.sum(y_axis*v_observed, axis=-1)
            pol_ang =  (initial_phase + np.arctan2(sinalpha, cosalpha)) % (2*np.pi)
        elif self.config.polarisation_modulation == "continuous_HWP":
            t_steps = self.get_t_steps()
            initial_phase = np.radians(self.config.HWP_phase_ini)
            HWP_angular_speed = 2*np.pi * self.config.HWP_rpm / 60.0         # RPM to radians/sec
            pol_angle = (initial_phase + t_steps*HWP_angular_speed) % (2*np.pi)
        elif self.config.polarisation_modulation == "stepped_HWP":
            t_steps = self.get_t_steps()
            initial_phase = np.radians(self.config.HWP_phase_ini)
            HWP_step = np.radians(self.config.HWP_step)                     # deg -> radians
            HWP_step_duration = self.config.HWP_step_duration
            pol_angle = (initial_phase + HWP_step*(t_step // HWP_step_duration)) % (2*np.pi)
        else:
            raise ValueError("polarisation modulation can be among the possible choices ['instrumental_scanning', 'continuous_HWP', 'stepped_HWP']. You have provided {}.".format(self.config.polarisation_modulation))
        return pol_ang

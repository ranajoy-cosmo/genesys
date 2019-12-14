import sys
import math
import numpy as np
import healpy as hp
import copy
from genesys.numerical.unit_conversion import Unit_Converter
import quaternion as qt
from genesys import Genesys_Class

uc = Unit_Converter()
t_year = uc.convert_unit(1.0, 'time', 'siderial year', 'second')

class Pointing(Genesys_Class):
    """
    Class for generating the pointing
    Left-handed Cartesian coordiante system
    The initial position of the satellite is assumed to be such that the anti-solar precession axis is the x-axis, 
    the spin axis is in the x-z plane and at a positive angle of alpha degrees to the precession axis, 
    the axis of revolution is the z-axis.
    """
    def __init__(self, pointing_params):
        self.copy_params(pointing_params)
        self.set_to_standard_units()
        self.set_rotation_axes()
        self.set_initial_axes_and_angles()
        self.set_initial_pointing_vector()

    def set_to_standard_units(self):
        self.params['alpha'] *= uc.conversion_factor('angle', 'degree', 'radian')
        self.params['beta'] *= uc.conversion_factor('angle', 'degree', 'radian')
        if 'pos' in self.params:
            self.params['pos'] *= uc.conversion_factor('angle', 'arcmin', 'radian')
        if 'HWP' in self.params:
            self.params['HWP']['rpm'] *= uc.conversion_factor('angular_speed', 'rpm', 'radians/sec')        # RPM to radians/sec

    def set_rotation_axes(self):
        self.axis_rev = np.array([0.0, 0.0, 1.0])       # z-axis, pointing upward
        self.axis_prec = np.array([1.0, 0.0, 0.0])      # x-axis, anti-solar axis
        self.axis_spin = np.array([np.cos(self.params['alpha']), 0.0, np.sin(self.params['alpha'])])        # at alpha degrees to axis_prec in the x-z plane

    def set_boresight_axes(self):
        boresight_opening_angle = self.params['alpha'] + self.params['beta']
        self.axis_boresight = np.array([np.cos(boresight_opening_angle), 0.0, np.sin(boresight_opening_angle)])
        self.axis_boresight_rot_x = np.array([np.cos(boresight_opening_angle - np.pi/2.0), 0.0, np.sin(boresight_opening_angle - np.pi/2.0)])

    def set_detector_axis(self):
        y_axis = np.array([0.0, 1.0, 0.0])

        q_rot_x = self.gen_rotation_quat(self.params['pos'][0], self.axis_boresight_rot_x)
        q_rot_y = self.gen_rotation_quat(self.params['pos'][1], -1.0*y_axis)

        q_rot_total = q_rot_x * q_rot_y

        self.detector_axis = self.rotate_vector(q_rot_total, self.axis_boresight)

    def set_detector_axis_offset(self, beam_integration_offset=[0.0,0.0]):
        # beam_integration_offset is given in arcmin
        y_axis = np.array([0.0, 1.0, 0.0])

        pos_offset = uc.convert_unit(self.params['pos'], 'angle', 'arcmin', 'radian') + uc.convert_unit(self.params['offset'], 'angle', 'arcsec', 'radian') + uc.convert_unit(beam_integration_offset, 'angle', 'arcmin', 'radian')

        q_rot_x = self.gen_rotation_quat(pos_offset, self.axis_boresight_rot_x)
        q_rot_y = self.gen_rotation_quat(pos_offset, -1.0*y_axis)

        q_rot_total = q_rot_x * q_rot_y

        self.detector_axis_offset = self.rotate_vector(q_rot_total, self.axis_boresight)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_t_steps(self, segment, segment_length, sampling_rate, n_samples):
        delta_t = 1.0 / sampling_rate
        t_start = (segment - 1) * segment_length
        t_stop = segment  * segment_length
        t_steps = t_start + delta_t * np.arange(n_samples)
        return t_steps

    def get_pointing_and_pol_angle(self, segment, segment_length, sampling_rate, n_samples, coordinate_system):
        """
        This is called for each new segment
        The simulation times, quaternions are set accordingly for the particular segment
        """
        t_steps = self.get_t_steps(segment, segment_length, sampling_rate, n_samples)
        self.total_rotation_quaternion = self.generate_quaternion(t_steps, coordinate_system)
        sky_pointing = self.rotate_vector(self.total_rotation_quaternion, self.pointing_vector_initial)
        pol_angle = self.get_pol_ang(t_steps, sky_pointing)
        theta, phi = hp.vec2ang(sky_pointing)
        return theta, phi, pol_angle

    def generate_quaternion(self, t_steps, coordinate_system):
        w_rev = 2 * np.pi / t_year
        w_prec = 2 * np.pi / self.params['scan_strategy']['t_precession']
        w_spin = 2 * np.pi / self.params['scan_strategy']['t_spin']
        total_rotation_quaternion = self.gen_rotation_quat(w_rev*t_steps, self.axis_rev) * self.gen_rotation_quat(w_prec*t_steps, self.axis_prec) * self.gen_rotation_quat(w_spin*t_steps, self.axis_spin)
        if coordinate_system == 'galactic':
            total_rotation_quaternion = self.quaternion_coordinate_transformation('E','G') * total_rotation_quaternion
        return total_rotation_quaternion

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    # Generating the polarisation angle
    def get_pol_ang(self, t_steps, v_observed):
        if self.params['pol_modulation'] == 'scan':
            initial_phase = np.radians(self.params['pol_phase_ini'])
            y_axis = self.rotate_vector(self.total_rotation_quaternion, np.array([0.0,1.0,0.0]))
            sinalpha = y_axis[...,0] * v_observed[...,1] - y_axis[...,1] * v_observed[...,0]
            cosalpha = y_axis[...,2] - v_observed[...,2] * np.sum(y_axis*v_observed, axis=-1)
            pol_ang =  (initial_phase + np.arctan2(sinalpha, cosalpha)) % (2*np.pi)
        else:
            HWP_angular_speed = 2*np.pi * self.params['HWP']['rpm'] / 60.0         # RPM to radians/sec
            pol_ang = (t_steps*HWP_angular_speed) % (2*np.pi)
        return pol_ang

    def get_hit_pix(self, theta, phi, nside):
        hit_pix = hp.ang2pix(nside, theta, phi)
        return hit_pix

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def gen_rotation_quat(self, thetas, axis_vec):
        """
        Initialising an Nx4 dimensional array of N quaternions of format (q_0, (q_1, q_2, q_3))

        Each row is an individual quaternion given by
        [cos(theta/2), sin(theta/2)*(x, y, z)]
        and, axis_vec=(x, y, z) is the axis of rotation.
        The axis_vec is an unit vector, and hence the quaternion will be an unit vector.
        Theta is in radians.
        The last two are considered to improve performance and save time on parsing additional lines.
        theta input shape : (N,) or a scalar     
        axis_vec input shape : (N,3) or (3,)
        output shape : (N,4) or (4,)
        """
        return qt.as_quat_array(np.insert(np.sin(thetas/2)[...,None]*axis_vec, 0, np.cos(thetas/2), axis=-1))

    def rotate_vector(self, rotation_quat, vec):
        """
        Rotate a given vector vec by the quaternion rotation_quat
        """
        quat_vec = qt.as_quat_array(np.insert(vec, 0, 0, axis=-1))
        return qt.as_float_array( rotation_quat * quat_vec * rotation_quat.conj() )[...,1:]

    def quaternion_coordinate_transformation(self, old_coord, new_coord, active=True):
        """
        Return the quaternion to rotate a vector in the old_coord frame to the new_coord frame
        """
        if active:
            transformation_euler = hp.rotator.get_coordconv_matrix([old_coord,new_coord])[0]
        else:
            transformation_euler = hp.rotator.get_coordconv_matrix([new_coord,old_coord])[0]
        return qt.from_rotation_matrix(transformation_euler)

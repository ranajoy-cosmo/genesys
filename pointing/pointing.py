import sys
import math
import numpy as np
import healpy as hp
import quaternion as qt
from genesys import Genesys_Class
from genesys.numerical.unit_conversion import Unit_Converter

uc = Unit_Converter()
t_year = uc.convert_unit(1.0, 'time', 'year', 'sec')
pi = np.pi

class Pointing(Genesys_Class):
    """
    Class for generating the pointing
    Left-handed Cartesian coordiante system
    By convention at time t=0, the precession(anti-solar) axis is along the x-axis
    The spin axis is in the x-z plane and at a positive angle of alpha degrees to the precession axis, 
    The axis of revolution is the z-axis.
    """
    def __init__(self, scan_strategy):
        self.params = {}
        self.params.update(scan_strategy)
        self.params['alpha'] *= uc.conversion_factor( 'angle', 'degree', 'radian')
        self.params['beta'] *= uc.conversion_factor('angle', 'degree', 'radian')
        self.set_rotation_axes()
        self.set_boresight_axis()

    def set_rotation_axes(self):
        self.axis_orbit = np.array([0.0, 0.0, 1.0])       # z-axis, pointing upward
        self.axis_prec = np.array([1.0, 0.0, 0.0])      # x-axis, anti-solar axis
        self.axis_spin = np.array([np.cos(self.params['alpha']), 0.0, np.sin(self.params['alpha'])])        # at alpha degrees to axis_prec in the x-z plane

    def set_boresight_axis(self):
        """
        The boresight axis is common for all detectors on the instrument
        It passes through the centre of the focal plane
        """
        boresight_opening_angle = self.params['alpha'] + self.params['beta']
        self.axis_boresight = np.array([np.cos(boresight_opening_angle), 0.0, np.sin(boresight_opening_angle)])

    def generate_obv_quaternion(self, t_steps, coordinate_system):
        """
        The observation quaternion is defined by 3 consecutive rotations
        1) Spin about axis_spin with angular speed w_spin
        2) Precession about axis_prec with angular speed w_prec
        3) Rotation about axis_orbit with angular speed w_orbit
        *4) Coordinate change to galactic, if specified
        """
        w_orbit = 2*pi / t_year
        w_prec = 2*pi / self.params['t_precession']
        w_spin = 2*pi / self.params['t_spin']
        total_rotation_quaternion = self.gen_rotation_quat(w_orbit*t_steps, self.axis_orbit) * self.gen_rotation_quat(w_prec*t_steps, self.axis_prec) * self.gen_rotation_quat(w_spin*t_steps, self.axis_spin)
        if coordinate_system == 'galactic':
            total_rotation_quaternion = self.quaternion_coordinate_transformation('E','G') * total_rotation_quaternion
        self.total_rotation_quaternion = total_rotation_quaternion

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_detector_axis(self, det_pos):
        """
        The pointing axis for each individual detector is generated here
        The axis for the polarisation axis is also generated here
        """
        _det_pos = np.array(det_pos) * uc.conversion_factor('angle', 'arcmin', 'radian')
        detector_opening_angle = self.params['alpha'] + self.params['beta'] + _det_pos[1]
        det_sight_vshifted = np.array([np.cos(detector_opening_angle), 0.0, np.sin(detector_opening_angle)])
        det_sight_orth_xz = np.array([np.cos(detector_opening_angle + pi/2.0), 0.0, np.sin(detector_opening_angle + pi/2.0)])
        q_rot_x = self.gen_rotation_quat(_det_pos[0], det_sight_orth_xz)
        axis_detector_sight = self.rotate_vector(q_rot_x, det_sight_vshifted)
        return axis_detector_sight

    def get_pointing_and_pol(self, axis_detector_sight):
        """
        Given the initial values of the boresight pointing the pointing for any subsequent set of t_steps can be determined. see next function for psi.
        """
        v_pointing = self.rotate_vector(self.total_rotation_quaternion, axis_detector_sight)
        theta, phi = hp.vec2ang(v_pointing)
        psi = self.get_polariser_angle(v_pointing)
        return theta, phi, psi

    def get_polariser_angle(self, v_pointing):
        """
        The psi is only for the intrinsic polarisation axis rotation 
        Along with the instrument motion. this value needs to be added to any initial orientation of the polarisation axis
        """
        y_axis = self.rotate_vector(self.total_rotation_quaternion, np.array([0.0,1.0,0.0]))
        sinalpha = y_axis[...,0] * v_pointing[...,1] - y_axis[...,1] * v_pointing[...,0]
        cosalpha = y_axis[...,2] - v_pointing[...,2] * np.sum(y_axis*v_pointing, axis=-1)
        psi =  np.arctan2(sinalpha, cosalpha)
        return psi

    def get_hwp_angle(self, t_steps, hwp_spin_rate):
        _hwp_spin_rate = uc.convert_unit(hwp_spin_rate, 'angular_velocity', 'rpm', 'radians/sec')
        hwp_psi = (_hwp_spin_rate * t_steps) % (2*np.pi)
        return hwp_psi

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # The Quaternion algebra section
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def gen_rotation_quat(self, thetas, axis_vec):
        """
        Initialising an Nx4 dimensional array of N quaternions of format (q_0, q_1, q_2, q_3)
        Each row is an individual quaternion given by
        [cos(theta/2), sin(theta/2)*(x, y, z)]
        and, axis_vec=(x, y, z) is the axis of rotation.
        The axis_vec is an unit vector, and hence the quaternion will be an unit vector.
        Theta is in radians.
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
        return qt.as_float_array(rotation_quat * quat_vec * rotation_quat.conj())[...,1:]

    def quaternion_coordinate_transformation(self, old_coord, new_coord, active=True):
        """
        Return the quaternion to rotate a vector in the old_coord frame to the new_coord frame
        """
        if active:
            transformation_euler = hp.rotator.get_coordconv_matrix([old_coord,new_coord])[0]
        else:
            transformation_euler = hp.rotator.get_coordconv_matrix([new_coord,old_coord])[0]
        return qt.from_rotation_matrix(transformation_euler)
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

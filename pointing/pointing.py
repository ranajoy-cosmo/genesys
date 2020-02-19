import sys
import math
import numpy as np
import healpy as hp
import copy
from genesys.numerical.unit_conversion import Unit_Converter
import quaternion as qt
from genesys import Genesys_Class

uc = Unit_Converter()
t_year = uc.convert_unit(1.0, 'time', 'siderial year', 'sec')
pi = np.pi

class Pointing(Genesys_Class):
    """
    
    CLASS FOR GENERATING THE POINTING
    LEFT-HANDED CARTESIAN COORDIANTE SYSTEM
    BY CONVENTION AT TIME t=0, THE PRECESSION(ANTI-SOLAR) AXIS IS ALONG THE x-axis
    THE SPIN AXIS IS IN THE X-Z PLANE AND AT A POSITIVE ANGLE OF alpha DEGREES TO THE PRECESSION AXIS, 
    THE AXIS OF REVOLUTION IS THE Z-AXIS.
    """
    def __init__(self, pointing_params):
        self.params = {}
        self.params.update(pointing_params)
        self.set_to_standard_units()
        self.set_rotation_axes()
        self.set_boresight_axis()
        self.set_detector_axis()

    def set_to_standard_units(self):
        self.params['alpha'] *= uc.conversion_factor('angle', 'degree', 'radian')
        self.params['beta'] *= uc.conversion_factor('angle', 'degree', 'radian')
        self.params['pos'] = np.array(self.params['pos'])
        self.params['pos'] *= uc.conversion_factor('angle', 'arcmin', 'radian')
        if 'HWP_spin_rate' in self.params:
            self.params['HWP_spin_rate'] *= uc.conversion_factor('angular_velocity', 'rpm', 'radians/sec')        # RPM to radians/sec

    def set_rotation_axes(self):
        self.axis_orbit = np.array([0.0, 0.0, 1.0])       # z-axis, pointing upward
        self.axis_prec = np.array([1.0, 0.0, 0.0])      # x-axis, anti-solar axis
        self.axis_spin = np.array([np.cos(self.params['alpha']), 0.0, np.sin(self.params['alpha'])])        # at alpha degrees to axis_prec in the x-z plane

    def set_boresight_axis(self):
        """
        THE BORESIGHT AXIS IS COMMON FOR ALL DETECTORS ON THE INSTRUMENT
        IT PASSES THROUGH THE CENTRE OF THE FOCAL PLANE
        """
        boresight_opening_angle = self.params['alpha'] + self.params['beta']
        self.axis_boresight = np.array([np.cos(boresight_opening_angle), 0.0, np.sin(boresight_opening_angle)])

    def set_detector_axis(self):
        """
        THE POINTING AXIS FOR EACH INDIVIDUAL DETECTOR IS GENERATED HERE
        THE AXIS FOR THE POLARISATION AXIS IS ALSO GENERATED HERE
        """
        detector_opening_angle = self.params['alpha'] + self.params['beta'] + self.params['pos'][1]
        det_sight_vshifted = np.array([np.cos(detector_opening_angle), 0.0, np.sin(detector_opening_angle)])
        det_sight_orth_xz = np.array([np.cos(detector_opening_angle + pi/2.0), 0.0, np.sin(detector_opening_angle + pi/2.0)])
        q_rot_x = self.gen_rotation_quat(self.params['pos'][0], det_sight_orth_xz)
        self.axis_detector_sight = self.rotate_vector(q_rot_x, det_sight_vshifted)
        # TO DO: POLARISAITION AXIS

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def generate_obv_quaternion(self, t_steps, coordinate_system):
        """
        THE OBSERVATION QUATERNION IS DEFINED BY 3 CONSECUTIVE ROTATIONS
        1) SPIN ABOUT AXIS_SPIN WITH ANGULAR SPEED W_SPIN
        2) PRECESSION ABOUT AXIS_PREC WITH ANGULAR SPEED W_PREC
        3) ROTATION ABOUT AXIS_ORBIT WITH ANGULAR SPEED W_ORBIT
        *4) COORDINATE CHANGE TO GALACTIC, IF SPECIFIED
        """
        # TO DO: Define inter-unit conversion factor method in unit_conversion
        w_orbit = 2*pi / t_year
        w_prec = 2*pi / self.params['t_precession']
        w_spin = 2*pi / self.params['t_spin']
        total_rotation_quaternion = self.gen_rotation_quat(w_orbit*t_steps, self.axis_orbit) * self.gen_rotation_quat(w_prec*t_steps, self.axis_prec) * self.gen_rotation_quat(w_spin*t_steps, self.axis_spin)
        if coordinate_system == 'galactic':
            total_rotation_quaternion = self.quaternion_coordinate_transformation('E','G') * total_rotation_quaternion
        return total_rotation_quaternion

    def get_pointing_and_pol_angle(self, t_steps, coordinate_system):
        """
        GIVEN THE INITIAL VALUES OF THE BORESIGHT POINTING THE POINTING FOR ANY SUBSEQUENT
        SET OF t_steps CAN BE DETERMINED. SEE NEXT FUNCTION FOR psi.
        """
        total_rotation_quaternion = self.generate_obv_quaternion(t_steps, coordinate_system)
        v_pointing = self.rotate_vector(total_rotation_quaternion, self.axis_detector_sight)
        theta, phi = hp.vec2ang(v_pointing)
        psi = self.get_polariser_angle(total_rotation_quaternion, v_pointing)
        return theta, phi, psi

    def get_polariser_angle(self, total_rotation_quaternion, v_pointing):
        """
        THE psi IS ONLY FOR THE INTRINSIC POLARISATION AXIS ROTATION 
        ALONG WITH THE INSTRUMENT MOTION. THIS VALUE NEEDS TO BE ADDED TO ANY INITIAL ORIENTATION
        OF THE POLARISATION AXIS
        """
        y_axis = self.rotate_vector(total_rotation_quaternion, np.array([0.0,1.0,0.0]))
        sinalpha = y_axis[...,0] * v_pointing[...,1] - y_axis[...,1] * v_pointing[...,0]
        cosalpha = y_axis[...,2] - v_pointing[...,2] * np.sum(y_axis*v_pointing, axis=-1)
        psi =  np.arctan2(sinalpha, cosalpha)
        return psi

    def get_hwp_angle(self, t_steps):
        hwp_psi = (self.params['HWP_spin_rate'] * t_steps) % (2*np.pi)
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

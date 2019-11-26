import numpy as np
import healpy as hp
import quaternion as qt
import genesys.numerical.quaternion as my_qt
import time

def gen_rotation_quat(thetas, axis_vec):
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

def rotate_vector(rotation_quat, vec):
    """
    Rotate a given vector vec by the quaternion rotation_quat
    """
    quat_vec = qt.as_quat_array(np.insert(vec, 0, 0, axis=-1))
    return qt.as_float_array( rotation_quat * quat_vec * rotation_quat.conj() )[...,1:]

def quaternion_coordinate_transformation(old_coord, new_coord, active=True):
    """
    Return the quaternion to rotate a vector in the old_coord frame to the new_coord frame
    """
    if active:
        transformation_euler = hp.rotator.get_coordconv_matrix([old_coord,new_coord])[0]
    else:
        transformation_euler = hp.rotator.get_coordconv_matrix([new_coord,old_coord])[0]
    return qt.from_rotation_matrix(transformation_euler)

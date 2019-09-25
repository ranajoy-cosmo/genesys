import numpy as np

def make_quaternion(thetas, vectors, degree=False):

    """
    Initialising an Nx4 dimensional array of N quaternions of format ndarray(w, x, y, z)

    Each row is an individual quaternion given by
    [cos(theta/2), sin(theta/2)*(x, y, z)]
    theta input shape : (N,)     
    vectors input shape : (N,3)
    output shape : (N,4)
    """

    if degree is True:
        thetas_rad = np.radians(thetas)
    else:
        thetas_rad=thetas
    return np.insert(np.sin(thetas_rad/2)[...,None]*vectors, 0, np.cos(thetas_rad/2), axis=-1)


def conjugate(q):
    """
    Finding the conjugate of an array of quaternions 

    Conjugate of (w, x, y, z) is (w, -x, -y, -z) 
    input shape : (N,4)
    output shape : (N,4)
    """

    return q*np.array([1, -1, -1, -1]) 


def norm(q):

    """
    Finding the norm of an array of quaternions 
    
    Norm of (w, x, y, z) is sqrt((w^2 + x^2 + y^2 + z^2)) 
    input Shape = (N,4)
    output shape = (N,1)
    """
    
    return np.sqrt(np.sum(q*q, axis=-1))[...,None]
    

def inverse(q):

    """
    Finding the inverse of an array of quaternions 

    Inverse of q is conjugate(q)/norm(q)^2 
    input Shape = (N,4)
    output shape = (N,4)
    """

    return conjugate(q)/(norm(q)**2) 
    

def normalise(q):

    """
    Normalising an array of quaternions 
    
    Norm of (w, x, y, z) is sqrt((w^2 + x^2 + y^2 + z^2)) 
    input shape : (N,4)
    output shape : (N,4)
    """

    return q/norm(q)


def vec2quat(vectors):

    """
    Construct an array of quaternions from an array of vectors with real part 0

    input shape : (N,3)
    output shape : (N,4)
    """

    return np.insert(vectors, 0, 0, axis=-1)


def multiply(q1, q2):
    
    """
    Multiplication of two quaternions a and b 

    input shape : (N,4), (N,4)
    output shape : (N,4)
    s -> scalar part
    v -> vector part
    a*b = (as*bs - av.bv, (as*bv + bs*av + avXbv))
    """
    
    q1s = q1[...,0]
    q2s = q2[...,0]
    q1v = q1[...,1:]
    q2v = q2[...,1:]

    qs = q1s*q2s - np.sum(q1v*q2v, axis=-1)
    qv = q1s[...,None]*q2v + q2s[...,None]*q1v + np.cross(q1v, q2v)

    return np.insert(qv, 0, qs, axis=-1) 


def transform(p, v, normalise=False):

    """
    Similarity Rotate a vector v by the normalised quaternion q

    v' = q*v*conj(q)
    input shape : (N,4), (N,3)
    output shape : (N,3)
    """

    if normalise is True:
        q = normalise(p)
    else:
        q = p

    return multiply(multiply(q, vec2quat(v)), conjugate(q))[...,1:]


def rotate_about_fixed_axis(thetas, axis, vectors, degree=False):
    
    """
    Rotate vectors about an axis by angles of theta 
    
    input shape : (N,), (3,), (N,3)
    output shape : (N,3)
    """

    return transform(make_quaternion(thetas, axis, degree), vectors)

def quaternion_from_euler(theta1, theta2, theta3, degree=False, mode='XYX'):

    """
    Create the quaternion corresponding to the Euler angles given
    """

    if mode is 'XYX':
        q = quaternion_XYX(theta1, theta2, theta3, degree)

    return q

def euler_rotations(theta1, theta2, theta3, vectors, degree=False, mode='XYX'):

    """
    Rotate a set of vectors by the 3 Euler angles given in the sequence mentioned
    """

    q = quaternion_from_euler(theta1, theta2, theta3, degree, mode)

    return transform(q, vectors)

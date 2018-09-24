#!/usr/bin/env python 

import numpy as np
from numpy import sin, cos
from numpy import hstack
from simulation.lib.geometry.conversions import deg2rad

"""
This module contains the conversion from Euler parameters to Quaternions.

The convention used here are those stated in the document 
Euler Angles, Quaternions and Transformation Matrices by D.M. Henderson 
for the NASA Shuttle Program.
"""


def quaternion_XYX(theta1, theta2, theta3, degree=False):

    if degree is True:
        theta1_rad = deg2rad(theta1)
        theta2_rad = deg2rad(theta2)
        theta3_rad = deg2rad(theta3)
    else:
        theta1_rad, theta2_rad, theta3_rad = theta1, theta2, theta3

    q1 = cos(theta2_rad/2)*cos((theta1_rad+theta3_rad)/2)
    q2 = cos(theta2_rad/2)*sin((theta1_rad+theta3_rad)/2)
    q3 = sin(theta2_rad/2)*cos((theta1_rad-theta3_rad)/2)
    q4 = sin(theta2_rad/2)*sin((theta1_rad-theta3_rad)/2)

    return hstack((q1[...,None], q2[...,None], q3[...,None], q4[...,None]))

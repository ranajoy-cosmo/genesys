import numpy as np
import healpy as hp
import os
import sys
from memory_profiler import profile


def get_dim(pol_type):
    """
    Calculates the dimension of a single pixel-pixel covariance matrix element
    and, the number of independent elements.
    """
    check_pol_type(pol_type)
    if pol_type == "IQU":
        dim = 3
    elif pol_type =="QU":
        dim = 2
    else:
        dim = 1

    ind_elements = int(dim*(dim+1)/2)

    return dim, ind_elements

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def flatten_block_matrix(block_matrix, pol_type):
    func_dict = {'IQU' : flatten_block_matrix_IQU, 'QU' : flatten_block_matrix_QU, 'I' : flatten_block_matrix_I}
    flat_matrix = func_dict[pol_type](block_matrix)
    return flat_matrix

def flatten_block_matrix_IQU(block_matrix):
    """
    In: 3x3 matrix.
    Out: 6 independent elements in a row.
    """
    num_elements = block_matrix.shape[0]
    flat_matrix = np.empty((num_elements, 6))

    flat_matrix[...,0] = block_matrix[..., 0, 0]
    flat_matrix[...,1] = block_matrix[..., 0, 0]
    flat_matrix[...,2] = block_matrix[..., 0, 2]
    flat_matrix[...,3] = block_matrix[..., 1, 1]
    flat_matrix[...,4] = block_matrix[..., 1, 2]
    flat_matrix[...,5] = block_matrix[..., 2, 2]

    return flat_matrix

def flatten_block_matrix_QU(block_matrix):
    """
    In: 2x2 matrix.
    Out: 3 independent elements in a row.
    """
    num_elements = block_matrix.shape[0]
    flat_matrix = np.empty((num_elements, 3))
    
    flat_matrix[..., 0] = block_matrix[..., 0, 0]
    flat_matrix[..., 1] = block_matrix[..., 0, 1]
    flat_matrix[..., 2] = block_matrix[..., 1, 1]

    return flat_matrix

def flatten_block_matrix_I(block_matrix):
    """
    In: 1 independent element for each pixel arranged in a 1-D array.
    Out: 1x1 matrix
    """
    num_elements = block_matrix.shape[0]
    flat_matrix = np.empty((num_elements, 1))

    flat_matrix[..., 0] = block_matrix[..., 0, 0]

    return flat_matrix

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def make_block_matrix(flat_matrix, pol_type):
    """
    Transforms a row of the independent elements in the covariance matrix to a square matrix.
    Works for multiple pixel-pixel elements
    """
    func_dict = {'IQU' : make_block_matrix_IQU, 'QU' : make_block_matrix_QU, 'I' : make_block_matrix_I}
    block_matrix = func_dict[pol_type](flat_matrix)
    return block_matrix

def make_block_matrix_IQU(flat_matrix):
    """
    In: 6 independent elements in a row.
    Out: 3x3 matrix.
    """
    num_blocks = flat_matrix.shape[0]

    block_matrix = np.empty((num_blocks, 3, 3))
    
    block_matrix[..., 0, 0] = flat_matrix[..., 0]
    block_matrix[..., 0, 1] = flat_matrix[..., 1]
    block_matrix[..., 0, 2] = flat_matrix[..., 2]
    block_matrix[..., 1, 0] = block_matrix[..., 0, 1] 
    block_matrix[..., 1, 1] = flat_matrix[..., 3]
    block_matrix[..., 1, 2] = flat_matrix[..., 4]
    block_matrix[..., 2, 0] = block_matrix[..., 0, 2] 
    block_matrix[..., 2, 1] = block_matrix[..., 1, 2] 
    block_matrix[..., 2, 2] = flat_matrix[..., 5]

    return block_matrix

def make_block_matrix_QU(flat_matrix):
    """
    In: 3 independent elements in a row.
    Out: 2x2 matrix.
    """
    num_blocks = flat_matrix.shape[0]

    block_matrix = np.empty((num_blocks, 2, 2))
    
    block_matrix[..., 0, 0] = flat_matrix[..., 0]
    block_matrix[..., 0, 1] = flat_matrix[..., 1]
    block_matrix[..., 1, 0] = block_matrix[..., 0, 1] 
    block_matrix[..., 1, 1] = flat_matrix[..., 2]

    return block_matrix

def make_block_matrix_I(flat_matrix):
    """
    In: 1 independent element for each pixel arranged in a 1-D array.
    Out: 1x1 matrix
    """
    num_blocks = flat_matrix.shape[0]

    block_matrix = np.empty((num_blocks, 1, 1))

    block_matrix[..., 0, 0] = flat_matrix[..., 0]

    return block_matrix

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def add_to_mm_matrices(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix, pol_type):
    """
    Generates the pixel-pixel inverse covariance matrix.
    This convention does not have a 1/2 multiplied in the data model. That factor is accounted for in the calibration.
    The inv_cov_matrix and b_matrix are passed by reference and are updated by the routine.
    """
    check_pol_type(pol_type)
    func_dict = {'IQU' : add_to_mm_matrices_IQU, 'QU' : add_to_mm_matrices_QU, 'I' : add_to_mm_matrices_I}
    func_dict[pol_type](hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix)
    
def add_to_mm_matrices_IQU(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix):
    n = np.bincount(hitpix, minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)

    hitmap += n

    inv_cov_matrix[..., 0] += n
    inv_cov_matrix[..., 1] += np.bincount(hitpix, weights=np.cos(2*pol), minlength=npix)
    inv_cov_matrix[..., 2] += np.bincount(hitpix, weights=np.sin(2*pol), minlength=npix)
    inv_cov_matrix[..., 3] += 0.5*(n + cos_4)
    inv_cov_matrix[..., 4] += 0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
    inv_cov_matrix[..., 5] += 0.5*(n - cos_4)

    b_matrix[..., 0] += np.bincount(hitpix, weights=ts, minlength=npix)
    b_matrix[..., 1] += np.bincount(hitpix, weights=ts*np.cos(2*pol), minlength=npix)
    b_matrix[..., 2] += np.bincount(hitpix, weights=ts*np.sin(2*pol), minlength=npix)

def add_to_mm_matrices_QU(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix):
    n = np.bincount(hitpix, minlength=npix)
    cos_4 = np.bincount(hitpix, weights=np.cos(4*pol), minlength=npix)

    hitmap += n

    inv_cov_matrix[..., 0] += 0.5*(n + cos_4)
    inv_cov_matrix[..., 1] += 0.5*np.bincount(hitpix, weights=np.sin(4*pol), minlength=npix)
    inv_cov_matrix[..., 2] += 0.5*(n - cos_4)

    b_matrix[..., 0] += np.bincount(hitpix, weights=ts*np.cos(2*pol), minlength=npix)
    b_matrix[..., 1] += np.bincount(hitpix, weights=ts*np.sin(2*pol), minlength=npix)

def add_to_mm_matrices_I(hitpix, pol, ts, inv_cov_matrix, b_matrix, hitmap, npix):
    n = np.bincount(hitpix, minlength=npix)

    hitmap += n

    inv_cov_matrix[..., 0] += n

    b_matrix[..., 0] += np.bincount(hitpix, weights=ts, minlength=npix)

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def get_covariance_matrix(inv_cov_matrix, hitmap, pol_type):
    dim, ind_elements = get_dim(pol_type)
    inv_cov_matrix_block = make_block_matrix(inv_cov_matrix, pol_type)
    mask_bad_pix(inv_cov_matrix_block, hitmap, pol_type)
    
    cov_matrix = np.linalg.inv(inv_cov_matrix_block)
    
    return cov_matrix

def get_sky_map(cov_matrix, b_matrix, hitmap, pol_type):
    sky_map = np.sum(cov_matrix*b_matrix[..., None], axis=1).T
    if pol_type == "IQU":
        sky_map[..., hitmap<3] = np.nan
    elif pol_type == "QU":
        sky_map[..., hitmap<2] = np.nan
    else:
        sky_map[..., hitmap<1] = np.nan
    return sky_map

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

def mask_bad_pix(inv_cov_matrix, hitmap, pol_type):
    if pol_type == "IQU":
        inv_cov_matrix[hitmap<3] = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    elif pol_type == "QU":
        inv_cov_matrix[hitmap<2] = np.array([[1.0, 0.0], [0.0, 1.0]])
    else:
        inv_cov_matrix[hitmap<1] = np.array([[1.0]])

def pointing_vec_to_hitpix(pointing_vec, nside):
    # Time-ordered array of pointing vectors to time-ordered array of hit pixels
    hitpix = hp.vec2pix(nside, pointing_vec[...,0], pointing_vec[...,1], pointing_vec[...,2])
    return hitpix

def write_maps(maps, map_type, field_names, recon_dir):

    hp.write_map(os.path.join(recon_dir, map_type+'.fits'), maps, column_names=field_names)

    if pol_type == "QU":
        hp.write_map(os.path.join(recon_dir, map_type+'.fits'), maps, column_names=['QQ', 'QU', 'UU'])
        #map_legends = {"QQ" : (0,0), "QU" : (0,1), "UU" : (1,1)}
    else:
        hp.write_map(os.path.join(recon_dir, map_type+'.fits'), maps, column_names=['TT', 'TQ', 'TU', 'QQ', 'QU', 'UU'])
        #map_legends = {"TT" : (0,0), "TQ" : (0,1), "TU" : (0,2), "QQ" : (1,1), "QU" : (1,2), "UU" : (2,2)}

    #for leg in map_legends.keys():
    #    hp.write_map(os.path.join(out_dir, "map_" + leg + ".fits"), maps[..., map_legends[leg][0], map_legends[leg][1]])

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
def check_pol_type(pol_type):
    assert pol_type in ['IQU', 'QU', 'I'], "pol_type has to be one of ['IQU','QU','I']. You have passes {}".format(pol_type)

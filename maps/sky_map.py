import numpy as np
import healpy as hp
from genesys import Genesys_Class
#  import genesys.noise.noise_utils as nu
import os

class Sky_Map(Genesys_Class):
    """
    THIS CLASS HAS USEFUL UTILITY FUNCTIONS FOR MANIPULATING HEALPix MAPS AND MASKS
    TWO MEMBER VARIABLES: sky_map, nside
    IF A sky_map IS MASKED, sky_map WILL CONTAIN THE MASKED MAP IN THE HEALPY FORMAT
    HEALPY MASK CONVENTION :
        0(False) -> PIXEL NOT SEEN
        1(True) -> PIXEL SEEN
    """
    def __init__(self, map_file_name=None, field=None, sky_map_np=None, other=None):
        """
        CONSTRUCTOR
        ORDER OF PREFERENCE:
            MAP_FILE_NAME -> READ_MAP_FROM_FILE(MAP_FILE_NAME, FIELD) : READ HEALPix MAP FROM FILE
            map_np -> from_np_array(map_np) : ASSIGN FROM A MAP EXISTING AS A NUMPY ARRAY
            other -> copy_attributes(other) : COPY CONSTRUCTOR
            None -> EMPTY OBJECT
        """
        if map_file_name is not None:
            self.read_map_from_file(map_file_name=map_file_name, field=field)
        elif sky_map_np is not None:
            self.from_np_array(sky_map_np)
        elif other is not None:
            self.copy_attributes(other=other)
        else:
            self.sky_map = None
            self.nside = None

    def read_map_from_file(self, map_file_name, field=None):
        """
        THE BASE DIRECTORY OF ALL MAPS IS ASSUMED TO BE global_paths['maps_dir']
        """
        if field == None:
            field = (0)
        self.sky_map = hp.read_map(self.get_path_to_map_file(map_file_name), field=field, verbose=False)
        self.nside = hp.get_nside(self.sky_map)

    def from_np_array(self, sky_map_np):
        """
        ACCEPTS A HEALPIX MAP IN NUMPY ARRAY FORMAT
        """
        self.sky_map = sky_map_np
        self.nside = hp.get_nside(self.sky_map)

    def write_map(self, map_file_name):
        """
        WRITES MAP TO FILE
        """
        hp.write_map(self.get_path_to_map_file(map_file_name), self.sky_map)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Path naming conventions
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def get_path_to_map_file(self, map_file_name):
        map_file_path = os.path.join(self.global_paths['maps_dir'], map_file_name)
        return map_file_path

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Mask routines
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def mask_map(self, binary_mask):
        """
        MASKS THE GIVEN SKY MAP WITH THE BINARY MASK.
        FOR EXAMPLE, PASSING A SINGLE BINARY MASK FOR A SET OF I,Q,U SKY MAPS, WILL APPLY THE SAME MASK TO ALL THE THREE.
        IN GENERAL, THE MASKS ARE ASSIGNED IN A CYCLIC MANNER. FOR EXAMPLE IF SKY_MAP HAS 4 MAPS AND BINARY_MASK HAS 2 MASKS, THEN:
        binary_mask[0] -> sky_map[0]
        binary_mask[1] -> sky_map[1]
        binary_mask[0] -> sky_map[2]
        binary_mask[1] -> sky_map[3]
        """
        # CHECKING THAT THE SKY MAP AND MASK HAVE THE SAME RESOLUTION
        assert self.nside == binary_mask.nside, "nside of sky map and mask does not match.\nnside of sky map : {}\nnside of mask : {}".format(self.nside, binary_mask.nside)

        self.sky_map = hp.ma(self.sky_map) 
        self.sky_map.mask = np.logical_not(binary_mask.sky_map)

    def get_mask_from_nan(self, labels=None):
        """
        Make a mask of the same dimension as the input sky mask.
        A NAN corresponds to an unseen pixel and will be set to 0 in the mask
        """
        if labels == None:
            labels = list(self.label_dict.keys())
        binary_mask = np.logical_not(self.sky_map == np.nan)
        boolean_selector = self.subset_label_selector(labels)
        return Sky_Map(sky_map_np=binary_mask[boolean_selector, :], map_labels=labels)
#
    #  def get_mask_from_cutoff(self, maximum=None, minimum=None, labels):
        #  if maximum == None:
            #  maximum = self.sky_map.max()
        #  if minimum == None:
            #  minimum = self.sky_map.min()
        #  binary_mask = np.ones(sky_mask.shape, dtype=np.int)
        #  binary_mask[sky_mask > maximum] = 0
        #  binary_mask[sky_mask < minimum] = 0
        #  return binary_mask
#
    #  def get_sky_fraction(binary_mask):
        #  """
        #  Returns the fraction of the sky that is valid.
        #  The function determines the dimension of the binary mask.
        #  It returns a scalar for a single mask. It returns an array of the same number of masks otherwise.
        #  """
        #  dim_mask = hp.maptype(binary_mask)
        #  n_pix = 12 * hp.get_nside(binary_mask)**2
        #  if dim_mask == 0:
            #  sky_fraction = np.sum(binary_mask.astype(np.float)) / n_pix
        #  else:
            #  sky_fraction = np.sum(binary_mask.astype(np.float), axis=1) / n_pix
        #  return sky_fraction
#
    #  def ud_grade_mask(mask, nside_out):
        #  """
        #  ud_grades the nside of the binary mask.
        #  Due to degrading, if a pixel has a value < 1, that is set to 0.
        #  """
        #  nside_in = hp.get_nside(mask)
        #  mask_new = hp.ud_grade(mask, nside_out=nside_out)
        #  if nside_out < nside_in:
            #  mask_new[mask_new < 1] = 0.0
        #  mask_new.dtype = np.int
#
        #  return mask_new

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    # Utilities
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def subset_label_selector(self, labels):
        boolean_selector = np.full(len(self.label_dict), False)
        for label in list(self.label_dict.keys()):
            if label in labels:
                boolean_selector[self.label_dict[label]] = True 
        return boolean_selector

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
#  def apodise_mask_gaussian(mask, fwhm, deg=True):
    #  # FWHM in degrees if deg==True, in arcmins otherwise
    #  if deg:
        #  mask_apodised = hp.smoothing(mask, fwhm=np.radians(fwhm))
    #  else:
        #  mask_apodised = hp.smoothing(mask, fwhm=np.radians(fwhm/60.0))
    #  return mask_apodised
#
#  def display_mask_statistics(mask):
    #  mask_dtype = mask.dtype
    #  mask_dim = hp.maptype(mask)
    #  nside = hp.get_nside(mask)
    #  sky_frac = get_sky_fraction(mask)
#
    #  print("#*#*#*")
    #  print("d-type\tdim\tnside\tsky-fraction")
    #  print("{}\t{}\t{}\t{}".format(mask_dtype, mask_dim, nside, sky_frac))
    #  print("#*#*#*\n")
#
#  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#
#  def fill_empty_pixels(sky_map, max_iter, fail_fill_value=0, pol=True, verbose=False):
    #  """
    #  Fill pixels with NAN with a fail_fill_value.
    #  """
    #  if np.sum(np.isnan(sky_map)) == 0:
        #  if verbose:
            #  prompt("There are no empty pixels")
        #  return
#
    #  nside = hp.get_nside(sky_map)
#
    #  if pol:
        #  dim = sky_map.shape[0]
        #  for i in xrange(max_iter):
            #  empty_pix = np.where(np.isnan(sky_map[0]))[0]
            #  theta, phi = hp.pix2ang(nside, empty_pix)
            #  neighbours = hp.get_all_neighbours(nside, theta, phi).T
            #  for j in range(dim):
                #  fill_values = np.nanmean(sky_map[j][neighbours], axis=-1)
                #  sky_map[j][empty_pix] = fill_values
            #  if np.sum(np.isnan(sky_map)) == 0:
                #  break
    #  else:
        #  for i in xrange(max_iter):
            #  empty_pix = np.where(np.isnan(sky_map))[0]
            #  theta, phi = hp.pix2ang(nside, empty_pix)
            #  neighbours = hp.get_all_neighbours(nside, theta, phi).T
            #  fill_values = np.nanmean(sky_map[neighbours], axis=-1)
            #  sky_map[empty_pix] = fill_values
            #  if np.sum(np.isnan(sky_map)) == 0:
                #  break
#
    #  num_empty_pix = np.sum(np.isnan(sky_map))
    #  if num_empty_pix:
        #  prompt("{} empty pixels remaining after {} iterations. Filling empty pixels with {}\n".format(num_empty_pix, max_iter, fail_fill_value))
        #  sky_map[np.isnan(sky_map)] = fail_fill_value
#
#  def deconvolve_map(map_in, fwhm_in=0.0, fwhm_out=0.0, lmax=None, binary_mask=None, pol=False, wiener=True, sky_prior=None):
    #
    #  if fwhm_in == fwhm_out:
        #  return map_in
#
    #  if lmax is None:
        #  lmax = 3*hp.get_nside(map_in) - 1
#
    #  if binary_mask is None:
        #  binary_mask = get_mask_from_map(map_in)
#
    #  f_sky = get_sky_fraction(binary_mask, pol)
#
    #  alm_in = su.estimate_alm(map_in, lmax, binary_mask, pol)
    #  alm_dec = su.deconvolve_alm(alm_in, fwhm_in=fwhm_in, fwhm_out=fwhm_out, f_sky=f_sky, pol=pol, wiener=True, sky_prior=sky_prior)
    #  map_dec = hp.alm2map(alm_dec, nside=hp.get_nside(map_in), pol=pol)
    #
    #  return map_dec

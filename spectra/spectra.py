import copy
import os
import numpy as np
import pandas as pd
import healpy as hp
import camb
from .. import Genesys_Class

class Spectra(Genesys_Class):
    """
    Class that handles I/O and manipulation of power spectra
    The spectra are handled as pandas DataFrame
    The multipole ell values are passed as the index and by default ranges between [0,lmax] for unbinned spectra
    The standard format in which specta are stored is [ell, TT, EE, BB, TE], and appropriate column names are applied anyway for additional columns, which may include TB, EB or anything as long as it's properly labelled
    The files are stored to disk as csv files
    The base directory of spectra files are by default global_paths.spectra_dir
    By default the Cl values and NOT the Dl values are written to file
    Conversion: l(l+1)Cl/2pi = Dl
    Member variables:
        spectra: Pandas DataFrame
    """

    # Constructor and I/O routines 
    """
    This block contains the following methods
        __init__
        read_spectra_from_file
        generate_spectra_from_camb_ini
        from_dict_spec
        write_spectra
    """
    def __init__(self, other=None, file_name=None, camb_ini_file_name=None, dict_spec=None, ells=None, lmax=None, columns=None):
        """
        The input parameters are:
            other: A Spectra object
            file_name: file from which the spectra in the prescribed format is to be read
            camb_ini: CAMB ini file used to generate new power spectra
            dict_spec: Spectra in the form of a dictionary
            ells, lmax and columns: Optional parameters
        The user must provide only one of the above parameters to get the desired results. Otherwise an empty object is created.
        """ 
        if other != None:
            self.copy_attributes(other=other)
        elif file_name != None:
            self.read_spectra_from_file(file_name=file_name, columns=columns,lmax=lmax)
        elif camb_ini_file_name != None:
            self.generate_spectra_from_camb_ini(camb_ini_file_name=camb_ini_file_name)
        elif dict_spec != None:
            self.from_dict_spec(dict_spec=dict_spec)
        else:
            self.spectra = pd.DataFrame(index=ells)
            self.spectra.index.name = 'ell'

    def read_spectra_from_file(self, file_name, lmax=None, columns=None):
        """
        Format for storing: csv
        The file name is provided with the appropriate extension
        """
        # TO DO: Check setting usecols for csv
        if columns != None:
            usecols = ["ell"] + columns
        else:
            usecols = None
        file_path = os.path.join(global_paths.spectra_dir, file_name)
        self.spectra = pd.read_csv(filepath_or_buffer=file_path, index_col="ell", usecols=usecols)
        if lmax:
            self.truncate(after=lmax)

    def generate_spectra_from_camb_ini(self, camb_ini_file_name):
        """
        Generates a new spectra from a CAMB ini file
        """
        camb_ini_file_path = os.path.join(self.global_paths['camb_params_dir'], camb_ini_file_name)
        params = camb.read_ini(camb_ini_file_path)
        results = camb.get_results(params)
        powers = results.get_cmb_power_spectra(params, CMB_unit='muK')
        self.from_dict_spec(dict_spec=powers)

    def from_dict_spec(self, dict_spec, ells=None, columns=None):
        """
        The spectra is loaded from a dictionary
        In this case it is assumed all the entries in the dictionary have common ells.
        If an entry for 'ell' is provided, it is used, otherwise it is assumed it is ordered from ell=0 to lmax
        """
        self.spectra = pd.DataFrame.from_dict(data=dict_spec)

    def write_spectra(self, relative_file_path):
        """
        Format for storing: csv
        The relative file path to global_paths.spectra_dir is provided as relative_file_path
        """
        spectra_file_path = os.path.join(self.global_paths.spectra_dir, relative_file_path)
        self.spectra.to_csv(path_or_buf=spectra_file_path, index_label="ell")

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    """
    This block contains methods for manipulating the ell values and are
        ells
        valid_ells_for_column
        invalid_ells_for_column
        lmax
        common_ells_with_nan
        add_empty_ell_rows
        remove_ells
        truncate
        sort_ells
        remove_ells_with_common_nan
        chop_trailing_ells_with_nans
        reindex_ell
    """

    def ells(self):
        """
        Return the ell index column for the dataframe<
        """
        return self.spectra.index.values

    def valid_ells_for_column(self, column):
        """
        Return the valid ells(not np.nan) for a particular column
        """
        return self.spectra.index[self.spectra[column].apply(np.isfinite)].values

    def invalid_ells_for_column(self, column):
        """
        Return all the ells for the column where it is np.nan
        """
        return self.spectra.index[self.spectra[column].apply(np.isnan)].values

    def lmax(self, column=None):
        if column == None:
            return self.ells[-1]
        else:
            return self.valid_ells_for_column(column)[-1]

    def common_ells_with_nan(self):
        common_ells = self.ells()
        for column in self.columns():
            common_ells = np.intersect1d(common_ells, self.invalid_ells_for_column(column))
        return common_ells

    def add_empty_ell_rows(self, ells):
        """
        Adds a row of ells filled with np.nan
        """
        columns = self.columns()
        col_size = len(columns)
        ell_size = len(ells)
        size = col_size * ell_size
        temp_df = pd.DataFrame(data=np.full(size, np.nan).reshape((ell_size,col_size)), index=ells, columns=columns)
        self.spectra = self.spectra.append(temp_df, verify_integrity=True)
        self.sort_ells()

    def remove_ells(self, ells):
        """
        Remove the given ells from the dataframe
        """
        self.spectra.drop(axis=0, index=ells, inplace=True)

    def truncate(self, before=None, after=None):
        """
        Chop off the spectra outside of before <= ell <= after
        The rows will be completely removed
        """
        self.spectra = self.spectra.truncate(before=before, after=after)

    def sort_ells(self):
        """
        Sort the ells in ascending order
        """
        self.spectra.sort_index(inplace=True)

    def remove_ells_with_common_nan(self):
        """
        Remove those ell indices which have np.nan in common with all columns
        """
        common_ells_with_nan = self.common_ells_with_nan()
        self.remove_ells(common_ells_with_nan)

    def chop_trailing_ells_with_nans(self):
        """
        Find the last ell value which is valid for all the columns
        Truncate the dataframe after the last valid ell
        """
        common_ells_with_nan = self.common_ells_with_nan()
        nan_groups = np.split(common_ells_with_nan, np.where(np.diff(common_ells_with_nan) != 1)[0]+1)
        self.remove_ells(nan_groups[-1])

    def reindex_ell(self, ells):
        """
        Force update the existing ells to the new set of ells.
        The existing values in columns are retained and new ells are filled with np.nan
        """
        self.spectra = self.spectra.reindex(index=ells, fill_value=np.nan, copy=True)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def columns(self):
        """
        Returns a list of the columns in the spectra dataframe
        """
        return self.spectra.columns.values

    def rename_columns(self, columns):
        """
        columns is a dictionary with {old_name: new_name}
        """
        self.spectra.rename(columns=columns, inplace=True)

    def reorder_columns(self, column_order):
        """
        Reorder the columns in self.spectra according to the new column_order
        """
        self.spectra = self.spectra[column_order]

    def common_columns(self, other):
        """
        Return list of common columns
        """
        return list(set(self.columns).intersection(other.columns))

    def uncommon_columns(self, other):
        """
        Return list of uncommon columns
        """
        return list(set(self.columns.symmetric_difference(other.columns)))

    def return_column_as_series(self, column, remove_nans=True):
        """
        Return the column in self as a pandas series
        """
        pass

    def add_empty_columns(self, columns):
        """
        Add a spectrum column with np.nan as empty value
        """
        for column in columns:
            if not column in self.spectra:
                self.spectra[column] = np.nan
            else:
                print("Column {} exists. Skipping. Try overwrite_empty_column()".format(column))

    def add_column(self, values, ells, column_name):
        """
        Add a new column labelled with the given ell
        Unspecified values will be filled by np.nan
        """
        new_column_df = pd.DataFrame(values, ells, columns=[column_name])
        self.spectra = self.spectra.join(new_column_df, how='outer')
        self.sort_ells()

    def remove_columns(self, columns):
        """
        Remove columns that are given in a list
        """
        self.spectra.drop(labels=columns, axis=1, inplace=True)

    def overwrite_with_empty_column(self, columns):
        """
        Overwrite a column with empty column containing np.nan
        """
        for column in columns:
            self.spectra[column] = np.nan

    def overwrite_column(self, values, ells, column_name):
        """
        Overwrite an existing column with the new column values
        """
        self.remove_columns(columns=[column_name])
        self.add_column(values, ells, column_name)

    def get_subset_columns(self, columns, return_new=True):
        """
        Keep only the subset of columns and remove the rest
        """
        drop_columns = set(self.columns()).symmetric_difference(columns) 
        if return_new:
            spectra_ret = Spectra(spectra_obj=self)
            spectra_ret.remove_columns(drop_columns)
            return spectra_ret
        else: 
            self.remove_columns(drop_columns)

    #  def join_spectra(self, other, columns_self=None, columns_other=None, return_new=True):
        #  """
        #  Join two spectra objects, i.e. their dataframe
        #  If columns_<> is None, all the columns are taken
        #  If return_new is True, a new object is returned keeping all the columns of self, otherwise they are added to self
        #  """
        #  if columns_self == None:
            #  columns_self = self.columns()
        #  if columns_other == None:
            #  columns_other = other.columns()
#
        #  if return_new:
            #  spectra_new = self.get_subset_columns(columns=columns_self)
            #  for
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def multiply(self, other, columns=None, return_new=False):
        """
        Multiply two spectra objects together
        If columns are given, only those columns are multiplied, otherwise the common columns are multiplied
        Only the ell values of the self object will be retained.
        ell values not present in other will result in np.nan in the result
        """
        if columns == None:
            columns = list(set(self.columns()).intersection(other.columns()))

        if return_new:
            ret_spectra = Spectra(spectra_obj=self)
            for column in columns:
                ret_spectra.spectra[column] *= other.spectra[column]
            return ret_spectra
        else:
            for column in columns:
                self.spectra[column] *= other.spectra[column]

    def divide(self, other, columns=None, return_new=False):
        """
        Multiply two spectra objects together
        If columns are given, only those columns are multiplied, otherwise the common columns are multiplied
        """
        if columns == None:
            columns = list(set(self.columns()).intersection(other.columns()))

        if return_new:
            ret_spectra = Spectra(spectra_obj=self)
            for column in columns:
                ret_spectra.spectra[column] /= other.spectra[column]
            return ret_spectra
        else:
            for column in columns:
                self.spectra[column] /= other.spectra[column]


    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def Cl_to_Dl(self, columns=None, return_new=False):
        """
        Transform from Cl to Dl
        Dl = l*(l+1)Cl/2*pi
        If return_obj is true, a new Spectra object is returned
        """
        ell = self.spectra.index.values
        if return_new:
            ret_spectra = Spectra(spectra_obj=self)
            ret_spectra.spectra = ret_spectra.spectra.multiply(ell*(ell+1) / (2*np.pi), axis='index')
            return ret_spectra
        else:
            self.spectra = self.spectra.multiply(ell*(ell+1) / (2*np.pi), axis='index')

    def Dl_to_Cl(self, columns=None, return_obj=False):
        """
        Transform from Dl to Cl
        Cl = Dl*2*pi/l*(l+1)
        If return_obj is true, a new Spectra object is returned
        """
        ell = self.spectra.index.values
        if return_obj:
            ret_spectra = Spectra(spectra_obj=self)
            ret_spectra.spectra = ret_spectra.spectra.multiply(2*np.pi / (ell*(ell+1)), axis='index')
            return ret_spectra
        else:
            self.spectra = self.spectra.multiply(2*np.pi / (ell*(ell+1)), axis='index')
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def generate_Bl_from_fwhm(self, fwhm, lmax, pol=True):
        """
        Get the window function assigned to self.spectra for a symmetric Gaussian beam
        fwhm is provided in arcmins and converted to radians internally
        """
        bl_np = hp.gauss_beam(fwhm=np.radians(fwhm/60.0), lmax=lmax, pol=pol)
        if pol:
            self.spectra = pd.DataFrame(bl_np, index=np.arange(0,lmax+1), columns=['TT','EE','BB','TE'])
        else:
            self.spectra = pd.DataFrame(bl_np, index=np.arange(0,lmax+1), columns=['TT'])

    def get_pixel_window_function(self, nside, lmax, pol=True):
        """
        Get the pixel window function assigned to self.spectra for the given nside
        """
        pixel_window_np = hp.pixwin(nside=nside, pol=pol)
        if pol:
            self.spectra['TT'] = pixel_window_np[0]
            self.spectra['EE'] = pixel_window_np[1]
            self.spectra['BB'] = pixel_window_np[1]
            self.spectra['TE'] = pixel_window_np[1]
        else:
            self.spectra['TT'] = pixel_window_np
        self.truncate(after=lmax)

    def convolve_spectra(self, columns, bl=None, fwhm=None, return_new=False):
        """
        bl takes precedence over fwhm
        bl is provided as a spectra object with columns corresponding to the power spectra
        fwhm is provided in arcmins and converted to radians internally
        bl is assumed to be an instance of class Spectra
        """
        if bl == None:
            bl = Spectra()
            bl.generate_Bl_from_fwhm(fwhm=fwhm, lmax=self.lmax())
        bl_squared = bl.multiply(bl)
        if return_new:
            ret_spectra = self.multiply(other=bl_squared, columns=columns, return_new=True)
            return ret_spectra
        else:
            self.multiply(other=bl_squared, columns=columns, return_new=False)

    def deconvolve_spectra(self, columns, bl=None, fwhm=None, return_new=False):
        """
        bl takes precedence over fwhm
        bl is provided as a spectra object with columns corresponding to the power spectra
        fwhm is provided in arcmins and converted to radians internally
        bl is assumed to be an instance of class Spectra
        """
        if bl == None:
            bl = Spectra()
            bl.generate_Bl_from_fwhm(fwhm=fwhm, lmax=self.lmax())
        bl_squared = bl.multiply(bl)
        if return_new:
            ret_spectra = self.divide(other=bl_squared, columns=columns, return_new=True)
            return ret_spectra
        else:
            self.divide(other=bl_squared, columns=columns, return_new=False)

    def pixel_window_convolve(self, nside, columns):
        """
        Convolve the spectra with the pixel window for the given nside
        """
        pixwin = Spectra()
        pixwin.get_pixel_window_function(nside=nside, lmax=2*nside, pol=True)
        self.convolve_spectra(columns=columns, bl=pixwin)

    def pixel_window_convolve(self, nside, columns):
        """
        Deconvolve the spectra with the pixel window for the given nside
        """
        pixwin = Spectra()
        pixwin.get_pixel_window_function(nside=nside, lmax=2*nside, pol=True)
        self.deconvolve_spectra(columns=columns, bl=pixwin)

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

    def bin_Dl(self, bin_edges=None, bin_centres_and_widths=None, columns=None, binned_spectra=None):
        """
        Bin the Dl in the given bins by simple average 
        bin_centres_and_widths take precedence over bin_edges
        The bin_edge is gievn as a list of tuples or 2-d array as [first_ell,last_ell]
            Bin centre is defined as (last_ell + first_ell)/2
            Bin width is defined as (last_ell - first_ell) + 1
        The bin_centres_and_widths is given as a list of tuples or 2-d array as [(bin_centre, bin_width)]
        The bin edges are then [bin_centre - int(bin_width/2), bin_centre + int(bin_width/2)]
        """
        if bin_centres_and_widths:
            bin_centres = np.array(bin_centres_and_widths)[...,0]
            bin_widths = np.array(bin_centres_and_widths)[...,1]
        else:
            first_ell = np.array(bin_edges)[...,0]
            last_ell = np.array(bin_edges)[...,1]
            bin_centres = (first_ell + last_ell) / 2
            bin_width = (last_ell - first_ell) + 1

        if binned_spectra == None:
            binned_spectra = Spectra()

        binned_column_temp = np.empty(bin_centres.size)
        for column in columns:
            for i in list(range(bin_centres.size)):
                binned_column_temp[i] = np.mean(self.spectra[column][first_ell,last_ell+1])
            binned_spectra.add_column(values=binned_column_temp, ells=bin_centres, column_name=column)
            binned_spectra.add_column(values=bin_widths, ells=bin_centres, column_name=column+"_width")

        return binned_spectra

    def bin_error_bar(error, bins, pol=False):
        n_bins = len(bins) - 1
        ell = np.arange(bins[-1] + 1)
        bin_width = np.empty(n_bins)
        if pol:
            binned_error = np.empty((3, n_bins))
        else:
            binned_error = np.empty(n_bins)

        for i in range(n_bins - 1):
            bin_low = bins[i]

    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

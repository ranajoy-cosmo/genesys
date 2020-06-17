<h1 align="center">Genesys</h3>
<p align="center">
  An end-to-end tool for simulating and analysing CMB space telescope data
</p>

## Features
* Generate **realistic pointing TODs** for a CMB space telescope having a rotation and precession, and revolution about the Sun.
* Generate **realistic signal TODs** containing Sky + Noise + Orbital dipole.
* **Simulate systematic** effects such as bandpass mismatch, pointing uncertainty....
* Instrument parameters containing scan strategy, channel parameters and individual detector parameters in easy to modifu yaml files.
* A **binning map-making algorithm** for estimation of I,Q,U Stokes parameters.
* A self-consistent and efficient I/O library built on HDF5.
* A Huffman compression system providing files in the Planck LFI format. 
* Highly parallelisable.

## Powered by
* numpy
* healpy
* scipy
* astropy
* pandas
* matplotlib
* quaternion
* camb
* h5py
* mpi4py
* ruamel.yaml

## Pipeline and data-flow

## Module description
* **tod**: The module for simulating the time-ordered instrument data.
* **map_maker**: Binning the TOD and relevant matrix manipulation routines.
* **pointing**: Generation of realistic instrument pointing and simple rotating HWP simulation.
* **gen_io**: Data I/O, Huffman compression, and data distribution.
* **instruments**: Parameter loading and manipulation methods for the instrument, channels and detector. Contains the configuraton for the instruments in respective sub-folders.
* **maps**: Loading, generation and manipulation of Healpix sky maps.
* **noise**: Generation of TOD noise and spectral analysis methods.
* **numerical**: Useful numerical modules including unit conversion, an alternate quaternion library, spectral analysis library.
* **physics**: Electromagnetic radiation physics methods and conversion between observational units.
* **plotting**: Set of plotting methods.
* **spectra**: Power spectum generation and manipulation methods.

## General setup

Clone the repo on to your local machine using

```
git clone https://github.com/ranajoy-cosmo/genesys.git
```

Go to the genesys parent directory with `cd genesys` and in the `__init__.py` file edit the `storage_dir = "/your/path/for/storing/data"`. This is the parent directory where all the big output files will be stored. Make sure you've sufficient disk space. Ideal for supercomputing clusters which provide TBs of space. 

To setup the storage directories and make soft links to them run `python setup_links.py link`. This will take care of making these directories and links. You will end up with:
* `genesys/output -> storage_dir/output`
* `genesys/data -> storage_dir/data`
* `genesys/maps/map_files -> storage_dir/maps`
* `genesys/spectra/spectra_files -> storage/spectra`

## Running the TOD simulation

This require setting up the instrument params.

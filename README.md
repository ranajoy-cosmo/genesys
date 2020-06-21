<h1 align="center">Genesys</h3>
<p align="center">
  An end-to-end tool for simulating and analysing CMB space telescope data
</p>

## Features
* Generate **realistic pointing TODs** for a CMB space telescope having a rotation about, a precession, and revolution about the Sun.
* Generate **realistic signal TODs** containing Sky + Noise + Orbital dipole + Systematics.
* Instrument parameters containing scan strategy, channel parameters and individual detector parameters in easy to modify yaml files.
* A **binning map-making algorithm** for estimation of I,Q,U Stokes parameters.
* A self-consistent and efficient I/O library built on **HDF5**.
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

This requires setting up a bit:
1. Populate the instrument and channel config files. A template is provided in `instruments/instrument_config/`
2. Populate the config file for the tod simulation run. A template is provided in `tod/config_files/default.yaml`

Relevant options for the parameters, units and descriptions are all proivided in the yaml files.

The sim_tod.py can be run in two modes: serial, and parallel.

For the serial run, do

```
python sim_tod.py <CONFIG_FILE> serial -v <VERBOSITY>
```

For the MPI run, do

```
mpirun -n <NUM_PROC> python sim_tod.py <CONFIG_FILE> mpi -v <VERBOSITY>
```

where, `CONFIG_FILE` is relative to the present directory, `VERBOSITY` is between 0 & 2, and `NUM_PROC` is the number of MPI processes.

## Data Model

The genesys data model is additive, adding in sky signal, noise, orbital dipole and any custom signal consecutively. The different components in the data model can be switched on/off by the corresponding parameters in the config file.

# Signal

The sky scan signal is given by 

![alt text](https://github.com/ranajoy-cosmo/genesys/blob/master/support_images/data_model.pdf)

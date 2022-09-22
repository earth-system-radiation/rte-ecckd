# rte-ecckd

This project aims to make the [ECCKD gas optics tool](https://doi.org/10.1029/2022MS003033) conform
to the radiative transfer API as defined in [RTE-RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp).

### How to build

##### Dependencies
The first step is to compile and install the following dependencies:
- A fortran compiler (gfortran, ifort, etc.)
- make
- [netCDF](https://github.com/Unidata/netcdf-c)
- [netCDF-fortran](https://github.com/Unidata/netcdf-fortran/)
- [hdf5](https://github.com/HDFGroup/hdf5)
- The [rte and rrtmgp libraries](https://github.com/earth-system-radiation/rte-rrtmgp)

##### Compiling the rte_ecckd library
Currently a simple `Makefile` is provided in the base of this repository.  This will be replaced with an
actual build system in the future.  For now, the library can be compiled using make:
```bash
FC=<path to fortran compiler> FCFLAGS="-I<directory where rte and rrtmgp modules are installed> -g -O0" make
```

##### Compiling the test programs
Test programs that run the [CMIP6 RFMIP RAD-IRF 100 column test cases](https://doi.org/10.5194/gmd-9-3447-2016) are provided.
They can also be compiled using make:
```bash
FC=<path to fortran compiler> \
FCFLAGS="-I<directory where rte and rrtmgp modules are installed> -I<directory where netcdff module is installed> -g -O0" \
LDFLAGS="-L. -L<directory where netcdff library is installed> -L<directory where netcdf library is installed> -L<directory where hdf5 library is installed>" \
LIBDIR="<directory where rte and rrtmgp static libraries are installed" \
make test
```

### How to run the tests
Before you can run the tests, you must download the necessay input data. A simple bash script is provided:
```bash
bash download-data-files.sh
```
After the data is downloaded, the tests can be run as follows:
```bash
# Calculates longwave fluxes
ecckd_rfmip_lw multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc data/ecckd-1.2_lw_ckd-definition_climate_fsck-tol0.0161.nc

# Calculate shortwave fluxes
ecckd_rfmip_sw multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc data/ecckd-1.2_sw_ckd-definition_climate_wide-tol0.05.nc
```




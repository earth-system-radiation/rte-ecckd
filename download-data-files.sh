#!/bin/bash

# Download the RFMIP RAD-IRF column data.
wget ftp://ftp.ldeo.columbia.edu/pub/robertp/rte-rrtmgp/continuous-integration/multiple_input4MIPs_radiation_RFMIP_UColorado-RFMIP-1-2_none.nc

# Download the template output files that will get filled with new data.  You must do this
# before running the model.
wget ftp://ftp.ldeo.columbia.edu/pub/robertp/rte-rrtmgp/continuous-integration/rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
mv rld_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rld_Efx_RTE-ecckd_rad-irf_r1i1p1f1_gn.nc

wget ftp://ftp.ldeo.columbia.edu/pub/robertp/rte-rrtmgp/continuous-integration/rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
mv rlu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rlu_Efx_RTE-ecckd_rad-irf_r1i1p1f1_gn.nc

wget ftp://ftp.ldeo.columbia.edu/pub/robertp/rte-rrtmgp/continuous-integration/rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
mv rsd_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rsd_Efx_RTE-ecckd_rad-irf_r1i1p1f1_gn.nc

wget ftp://ftp.ldeo.columbia.edu/pub/robertp/rte-rrtmgp/continuous-integration/rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc
mv rsu_Efx_RTE-RRTMGP-181204_rad-irf_r1i1p1f1_gn.nc rsu_Efx_RTE-ecckd_rad-irf_r1i1p1f1_gn.nc

program ecckd_rfmip_lw
use, intrinsic :: iso_fortran_env, only: error_unit
use gas_optics_ecckd, only: ty_gas_optics_ecckd
use load_coefficients, only: load_and_init
use mo_fluxes, only: ty_fluxes_broadband
use mo_gas_concentrations, only: ty_gas_concs
use mo_optical_props, only: ty_optical_props_1scl
use mo_rte_kind, only: wp
use mo_rte_lw, only: rte_lw
use mo_source_functions, only: ty_source_func_lw
use rfmip_io, only: read_and_block_gases_ty, read_and_block_lw_bc, read_and_block_pt, &
                    read_size, unblock_and_write
use simple_netcdf, only: stop_on_err
use utils, only: determine_gas_names, parse_args
implicit none


character(len=132) :: rfmip_path, ecckd_path
character(len=132) :: flxdn_file, flxup_file
integer :: ncol, nlay, nbnd, nexp, nblocks, block_size, &
           forcing_index, physics_index, n_quad_angles
logical :: top_at_1
integer :: b, icol, ibnd
character(len=1) :: forcing_index_char, physics_index_char
character(len=32), dimension(6) :: names_in_kdist, names_in_rfmip
real(kind=wp), dimension(:,:,:), allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
real(kind=wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
real(kind=wp), dimension(:,:), allocatable :: sfc_emis, sfc_t ! block_size, nblocks (emissivity is spectrally constant)
real(kind=wp), dimension(:,:), allocatable :: sfc_emis_spec ! nbands, block_size (spectrally-resolved emissivity)
type(ty_gas_optics_ecckd) :: ecckd
type(ty_source_func_lw) :: source
type(ty_optical_props_1scl) :: optical_props
type(ty_fluxes_broadband) :: fluxes
type(ty_gas_concs), dimension(:), allocatable :: gas_conc_array


!Parse command line arguments.
call parse_args(rfmip_path, ecckd_path, forcing_index, physics_index)
block_size = 1
if (physics_index .eq. 2) then
  n_quad_angles = 3
else
  n_quad_angles = 1
endif

!How big is the problem? Does it fit into blocks of the size we've specified?
call read_size(rfmip_path, ncol, nlay, nexp)
if (mod(ncol*nexp, block_size) .ne. 0) then
  call stop_on_err("ecckd_rfmip_lw: number of columns doesn't fit evenly into blocks.")
endif
nblocks = (ncol*nexp)/block_size
write(error_unit, *) "Using ",  nblocks, " blocks of size ", block_size

!Define output files.
write(forcing_index_char, "(i1)") forcing_index
write(physics_index_char, "(i1)") physics_index
write(error_unit, *) "Using forcing index "//trim(forcing_index_char)// &
                     " and physics index "//trim(physics_index_char)
flxdn_file = "rld_Efx_RTE-ecckd_rad-irf_r1i1p" //  &
             trim(physics_index_char) // "f" // trim(forcing_index_char) // "_gn.nc"
flxup_file = "rlu_Efx_RTE-ecckd_rad-irf_r1i1p" // &
             trim(physics_index_char) // "f" // trim(forcing_index_char) // "_gn.nc"

!Identify the set of gases used in the calculation based on the forcing index
!A gas might have a different name in the k-distribution than in the files
!provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
call determine_gas_names(forcing_index, names_in_kdist, names_in_rfmip)
write(error_unit, *) "Calculation uses RFMIP gases: ", &
                     (trim(names_in_rfmip(b))//" ", b=1, size(names_in_rfmip))

!Read RFMIP input data.
call read_and_block_pt(rfmip_path, block_size, p_lay, p_lev, t_lay, t_lev)
call read_and_block_lw_bc(rfmip_path, block_size, sfc_emis, sfc_t)
call read_and_block_gases_ty(rfmip_path, block_size, names_in_kdist, &
                             names_in_rfmip, gas_conc_array)

!Read ecckd data.
call load_and_init(ecckd, trim(ecckd_path), gas_conc_array(1))
if (.not. ecckd%source_is_internal() ) then
  call stop_on_err("ecckd_rfmip_lw: k-distribution file isn't for longwave.")
endif
nbnd = ecckd%get_nband()

!Which way is up?
top_at_1 = p_lay(1, 1, 1) .lt. p_lay(1, nlay, 1)

!RRTMGP won't run with pressure less than its minimum. The top level in the RFMIP file
!is set to 10^-3 Pa. Here we pretend the layer is just a bit less deep.
!This introduces an error but shows input sanitizing.
if (top_at_1) then
  p_lev(:,1,:) = ecckd%get_press_min() + epsilon(ecckd%get_press_min())
else
  p_lev(:,nlay+1,:) = ecckd%get_press_min() + epsilon(ecckd%get_press_min())
endif

!Allocate space for output fluxes (accessed via pointers in ty_fluxes_broadband),
!gas optical properties, and source functions. The %alloc() routines carry along
!the spectral discretization from the k-distribution.
allocate(flux_up(block_size, nlay+1, nblocks), &
         flux_dn(block_size, nlay+1, nblocks))
allocate(sfc_emis_spec(nbnd, block_size))
call stop_on_err(source%alloc(block_size, nlay, ecckd))
call stop_on_err(optical_props%alloc_1scl(block_size, nlay, ecckd))

!Loop over blocks
!do b = 1, nblocks
do b = 1, 1700
  fluxes%flux_up => flux_up(:,:,b)
  fluxes%flux_dn => flux_dn(:,:,b)

  !Expand the spectrally-constant surface emissivity to a per-band emissivity for each column
  do icol = 1, block_size
    do ibnd = 1, nbnd
      sfc_emis_spec(ibnd,icol) = sfc_emis(icol,b)
    enddo
  enddo

  !Compute the optical properties of the atmosphere and the Planck source functions
  !from pressures, temperatures, and gas concentrations...
  call stop_on_err(ecckd%gas_optics(p_lay(:,:,b), &
                                    p_lev(:,:,b), &
                                    t_lay(:,:,b), &
                                    sfc_t(:,b), &
                                    gas_conc_array(b), &
                                    optical_props, &
                                    source, &
                                    tlev=t_lev(:,:,b)))

  !Compute the spectrally-resolved fluxes, providing reduced values via ty_fluxes_broadband.
  call stop_on_err(rte_lw(optical_props, &
                          top_at_1, &
                          source, &
                          sfc_emis_spec, &
                          fluxes, &
                          n_gauss_angles=n_quad_angles))
enddo

!Write out the output.
call unblock_and_write(trim(flxup_file), "rlu", flux_up)
call unblock_and_write(trim(flxdn_file), "rld", flux_dn)


end program ecckd_rfmip_lw

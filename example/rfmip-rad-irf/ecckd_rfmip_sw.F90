program ecckd_rfmip_sw
use, intrinsic :: iso_fortran_env, only: error_unit
use mo_rte_kind, only: wp
!use mo_rte_util_array, only: zero_array
use mo_optical_props, only: ty_optical_props_2str
use gas_optics_ecckd, only: ty_gas_optics_ecckd
use mo_gas_concentrations, only: ty_gas_concs
use mo_rte_sw, only: rte_sw
use mo_fluxes, only: ty_fluxes_broadband
use mo_load_coefficients,  only: load_and_init
use mo_rfmip_io, only: read_size, read_and_block_pt, read_and_block_gases_ty, &
                       unblock_and_write, read_and_block_sw_bc
use mo_simple_netcdf, only: stop_on_err
use utils, only: determine_gas_names, parse_args
implicit none


character(len=132) :: rfmip_path, ecckd_path
character(len=132) :: flxdn_file, flxup_file
integer :: ncol, nlay, nbnd, ngpt, nexp, nblocks, block_size, forcing_index, physics_index
logical :: top_at_1
integer :: b, icol, ibnd, igpt
character(len=1) :: forcing_index_char, physics_index_char
character(len=32), dimension(6) :: names_in_kdist, names_in_rfmip
real(wp), dimension(:,:,:), allocatable :: p_lay, p_lev, t_lay, t_lev ! block_size, nlay, nblocks
real(wp), dimension(:,:,:), target, allocatable :: flux_up, flux_dn
real(wp), dimension(:,:), allocatable :: surface_albedo, total_solar_irradiance, &
                                         solar_zenith_angle
real(wp), dimension(:,:), allocatable :: sfc_alb_spec
type(ty_gas_optics_ecckd) :: ecckd
type(ty_optical_props_2str) :: optical_props
type(ty_fluxes_broadband) :: fluxes
real(wp), dimension(:,:), allocatable :: toa_flux ! block_size, ngpt
real(wp), dimension(:), allocatable :: def_tsi, mu0    ! block_size
logical, dimension(:,:), allocatable :: usecol ! block_size, nblocks
type(ty_gas_concs), dimension(:), allocatable :: gas_conc_array
real(wp), parameter :: deg_to_rad = acos(-1._wp)/180._wp


!Parse command line arguments.
call parse_args(rfmip_path, ecckd_path, forcing_index, physics_index)
block_size = 1

!How big is the problem? Does it fit into blocks of the size we've specified?
call read_size(rfmip_path, ncol, nlay, nexp)
if (mod(ncol*nexp, block_size) .ne. 0) then
  call stop_on_err("rrtmgp_rfmip_sw: number of columns doesn't fit evenly into blocks.")
endif
nblocks = (ncol*nexp)/block_size
write(error_unit, *) "Using ",  nblocks, " blocks of size ", block_size

!Define output files.
write(forcing_index_char, "(i1)") forcing_index
write(physics_index_char, "(i1)") physics_index
write(error_unit, *) "Using forcing index "//trim(forcing_index_char)// &
                     " and physics index "//trim(physics_index_char)
flxdn_file = "rsd_Efx_RTE-ecckd_rad-irf_r1i1p1f"//trim(forcing_index_char)//"_gn.nc"
flxup_file = "rsu_Efx_RTE-ecckd_rad-irf_r1i1p1f"//trim(forcing_index_char)//"_gn.nc"

!Identify the set of gases used in the calculation based on the forcing index
!A gas might have a different name in the k-distribution than in the files
!provided by RFMIP (e.g. 'co2' and 'carbon_dioxide')
call determine_gas_names(forcing_index, names_in_kdist, names_in_rfmip)
write(error_unit, *) "Calculation uses RFMIP gases: ", &
                     (trim(names_in_rfmip(b))//" ", b=1, size(names_in_rfmip))

!Read RFMIP input data.
call read_and_block_pt(rfmip_path, block_size, p_lay, p_lev, t_lay, t_lev)
call read_and_block_sw_bc(rfmip_path, block_size, surface_albedo, &
                          total_solar_irradiance, solar_zenith_angle)
call read_and_block_gases_ty(rfmip_path, block_size, names_in_kdist, &
                             names_in_rfmip, gas_conc_array)

!Read ecckd data.
call load_and_init(ecckd, trim(ecckd_path), gas_conc_array(1))
if (.not. ecckd%source_is_external()) then
  call stop_on_err("ecckd_rfmip_sw: k-distribution file isn't for shortwave.")
endif
nbnd = ecckd%get_nband()
ngpt = ecckd%get_ngpt()

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
allocate(mu0(block_size), sfc_alb_spec(nbnd, block_size))
allocate(toa_flux(block_size, ngpt), &
         def_tsi(block_size), usecol(block_size, nblocks))
call stop_on_err(optical_props%alloc_2str(block_size, nlay, ecckd))

!RTE will fail if passed solar zenith angles greater than 90 degree. We replace any with
!nighttime columns with a default solar zenith angle. We'll mask these out later, of
!course, but this gives us more work and so a better measure of timing.
do b = 1, nblocks
  usecol(1:block_size,b) = solar_zenith_angle(1:block_size,b) .lt. 90._wp - 2._wp*spacing(90._wp)
enddo

!Loop over blocks
!do b = 1, nblocks
do b = 1, 1700
  fluxes%flux_up => flux_up(:,:,b)
  fluxes%flux_dn => flux_dn(:,:,b)

  !Compute the optical properties of the atmosphere and the Planck source functions
  !from pressures, temperatures, and gas concentrations...
  call stop_on_err(ecckd%gas_optics(p_lay(:,:,b), &
                                    p_lev(:,:,b), &
                                    t_lay(:,:,b), &
                                    gas_conc_array(b), &
                                    optical_props, &
                                    toa_flux))
  !Boundary conditions
  !What's the total solar irradiance assumed by RRTMGP?
  def_tsi(1:block_size) = sum(toa_flux, dim=2)

  !Normalize incoming solar flux to match RFMIP specification
  do igpt = 1, ngpt
    do icol = 1, block_size
      toa_flux(icol,igpt) = toa_flux(icol,igpt)*total_solar_irradiance(icol,b)/def_tsi(icol)
    enddo
  enddo

  !Expand the spectrally-constant surface albedo to a per-band albedo for each column
  do icol = 1, block_size
    do ibnd = 1, nbnd
      sfc_alb_spec(ibnd,icol) = surface_albedo(icol,b)
    enddo
  enddo

  !Cosine of the solar zenith angle
  do icol = 1, block_size
    mu0(icol) = merge(cos(solar_zenith_angle(icol,b)*deg_to_rad), 1._wp, usecol(icol,b))
  enddo

  !Compute the spectrally-resolved fluxes, providing reduced values via ty_fluxes_broadband.
  call stop_on_err(rte_sw(optical_props, &
                          top_at_1, &
                          mu0, &
                          toa_flux, &
                          sfc_alb_spec, &
                          sfc_alb_spec, &
                          fluxes))
  !Zero out fluxes for which the original solar zenith angle is > 90 degrees.
  do icol = 1, block_size
    if (.not. usecol(icol,b)) then
      flux_up(icol,:,b)  = 0._wp
      flux_dn(icol,:,b)  = 0._wp
    endif
  enddo
enddo

!Write out the output.
call unblock_and_write(trim(flxup_file), "rsu", flux_up)
call unblock_and_write(trim(flxdn_file), "rsd", flux_dn)


end program ecckd_rfmip_sw

! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! The gas optics class used by RRMTGP needs to be initialized with data stored in a netCDF file.
!    RRTMGP itself doesn't include methods for reading the data so we don't conflict with users'
!    local environment. This module provides a straight-forward implementation of reading the data
!    and calling gas_optics%load().
!
! -------------------------------------------------------------------------------------------------
module mo_load_coefficients
  !
  ! Modules for working with rte and rrtmgp
  !
  use mo_rte_kind,           only: wp, wl
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_optics_ecckd, only: AbsorptionTable, none_, linear, look_up_table, &
                                 relative_linear, ty_gas_optics_ecckd
  ! --------------------------------------------------
  use mo_simple_netcdf, only: read_field, read_char_vec, read_logical_vec, var_exists, &
                              get_dim_size, get_var_size
  use netcdf
  implicit none
  private
  public :: load_and_init


  interface load_and_init
    module procedure load_and_init_rrtmgp
    module procedure load_and_init_ecckd
  end interface


contains
  subroutine stop_on_err(msg)
    use iso_fortran_env, only : error_unit
    character(len=*), intent(in) :: msg


    if(msg /= "") then
      write(error_unit, *) msg
      error stop 1
    end if
  end subroutine
  !--------------------------------------------------------------------------------------------------------------------
  ! read optical coefficients from NetCDF file
  subroutine load_and_init_rrtmgp(kdist, filename, available_gases)
    class(ty_gas_optics_rrtmgp), intent(inout) :: kdist
    character(len=*),     intent(in   ) :: filename
    class(ty_gas_concs),  intent(in   ) :: available_gases ! Which gases does the host model have available?
    ! --------------------------------------------------
    !
    ! Variables that will be passed to gas_optics%load()
    !
    character(len=32), dimension(:), allocatable :: gas_names
    integer,  dimension(:,:,:),      allocatable :: key_species
    integer,  dimension(:,:  ),      allocatable :: band2gpt
    real(wp), dimension(:,:  ),      allocatable :: band_lims
    real(wp)                                     :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:      ),    allocatable :: press_ref
    real(wp), dimension(:      ),    allocatable :: temp_ref
    real(wp), dimension(:,:,:  ),    allocatable :: vmr_ref
    real(wp), dimension(:,:,:,:),    allocatable :: kmajor

    character(len=32), dimension(:),  allocatable :: gas_minor, identifier_minor
    character(len=32), dimension(:),  allocatable :: minor_gases_lower,               minor_gases_upper
    integer, dimension(:,:),          allocatable :: minor_limits_gpt_lower,          minor_limits_gpt_upper
    logical(wl), dimension(:),        allocatable :: minor_scales_with_density_lower, minor_scales_with_density_upper
    character(len=32), dimension(:),  allocatable :: scaling_gas_lower,               scaling_gas_upper
    logical(wl), dimension(:),        allocatable :: scale_by_complement_lower,       scale_by_complement_upper
    integer, dimension(:),            allocatable :: kminor_start_lower,              kminor_start_upper
    real(wp), dimension(:,:,:),       allocatable :: kminor_lower,                    kminor_upper

    real(wp), dimension(:,:,:  ), allocatable :: rayl_lower, rayl_upper
    real(wp), dimension(:      ), allocatable :: solar_quiet, solar_facular, solar_sunspot
    real(wp)                                  :: tsi_default, mg_default, sb_default
    real(wp), dimension(:,:    ), allocatable :: totplnk
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac
    real(wp), dimension(:,:)    , allocatable :: optimal_angle_fit

    ! -----------------
    !
    ! Book-keeping variables
    !
    integer :: ncid
    integer :: ntemps,          &
               npress,          &
               nabsorbers,      &
               nextabsorbers,   &
               nminorabsorbers, &
               nmixingfracs,    &
               nlayers,         &
               nbnds,           &
               ngpts,           &
               npairs,          &
               nminor_absorber_intervals_lower, &
               nminor_absorber_intervals_upper, &
               ncontributors_lower, &
               ncontributors_upper, &
               ninternalSourcetemps, &
               nfit_coeffs
    ! --------------------------------------------------
    !
    ! How big are the various arrays?
    !
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("load_and_init(): can't open file " // trim(fileName))
    ntemps            = get_dim_size(ncid,'temperature')
    npress            = get_dim_size(ncid,'pressure')
    nabsorbers        = get_dim_size(ncid,'absorber')
    nminorabsorbers   = get_dim_size(ncid,'minor_absorber')
    nextabsorbers     = get_dim_size(ncid,'absorber_ext')
    nmixingfracs      = get_dim_size(ncid,'mixing_fraction')
    nlayers           = get_dim_size(ncid,'atmos_layer')
    nbnds             = get_dim_size(ncid,'bnd')
    ngpts             = get_dim_size(ncid,'gpt')
    npairs            = get_dim_size(ncid,'pair')
    nminor_absorber_intervals_lower &
                      = get_dim_size(ncid,'minor_absorber_intervals_lower')
    nminor_absorber_intervals_upper  &
                      = get_dim_size(ncid,'minor_absorber_intervals_upper')
    ninternalSourcetemps &
                      = get_dim_size(ncid,'temperature_Planck')
    ncontributors_lower = get_dim_size(ncid,'contributors_lower')
    ncontributors_upper = get_dim_size(ncid,'contributors_upper')
    nfit_coeffs         = get_dim_size(ncid,'fit_coeffs') ! Will be 0 for SW

    ! -----------------
    !
    ! Read the many arrays
    !
    gas_names         = read_char_vec(ncid, 'gas_names', nabsorbers)
    key_species       = int(read_field(ncid, 'key_species',  2, nlayers, nbnds))
    band_lims         = read_field(ncid, 'bnd_limits_wavenumber', 2, nbnds)
    band2gpt          = int(read_field(ncid, 'bnd_limits_gpt', 2, nbnds))
    press_ref         = read_field(ncid, 'press_ref', npress)
    temp_ref          = read_field(ncid, 'temp_ref',  ntemps)
    temp_ref_p        = read_field(ncid, 'absorption_coefficient_ref_P')
    temp_ref_t        = read_field(ncid, 'absorption_coefficient_ref_T')
    press_ref_trop    = read_field(ncid, 'press_ref_trop')
    kminor_lower      = read_field(ncid, 'kminor_lower', &
        ncontributors_lower, nmixingfracs, ntemps)
    kminor_upper      = read_field(ncid, 'kminor_upper', &
        ncontributors_upper, nmixingfracs, ntemps)
    gas_minor = read_char_vec(ncid, 'gas_minor', nminorabsorbers)
    identifier_minor = read_char_vec(ncid, 'identifier_minor', nminorabsorbers)
    minor_gases_lower = read_char_vec(ncid, 'minor_gases_lower', nminor_absorber_intervals_lower)
    minor_gases_upper = read_char_vec(ncid, 'minor_gases_upper', nminor_absorber_intervals_upper)
    minor_limits_gpt_lower &
                      = int(read_field(ncid, 'minor_limits_gpt_lower', npairs,nminor_absorber_intervals_lower))
    minor_limits_gpt_upper &
                      = int(read_field(ncid, 'minor_limits_gpt_upper', npairs,nminor_absorber_intervals_upper))
    minor_scales_with_density_lower &
                      = read_logical_vec(ncid, 'minor_scales_with_density_lower', nminor_absorber_intervals_lower)
    minor_scales_with_density_upper &
                      = read_logical_vec(ncid, 'minor_scales_with_density_upper', nminor_absorber_intervals_upper)
    scale_by_complement_lower &
                      = read_logical_vec(ncid, 'scale_by_complement_lower', nminor_absorber_intervals_lower)
    scale_by_complement_upper &
                      = read_logical_vec(ncid, 'scale_by_complement_upper', nminor_absorber_intervals_upper)
    scaling_gas_lower &
                      = read_char_vec(ncid, 'scaling_gas_lower', nminor_absorber_intervals_lower)
    scaling_gas_upper &
                      = read_char_vec(ncid, 'scaling_gas_upper', nminor_absorber_intervals_upper)
    kminor_start_lower &
                      = int(read_field(ncid, 'kminor_start_lower', nminor_absorber_intervals_lower))
    kminor_start_upper &
                      = int(read_field(ncid, 'kminor_start_upper', nminor_absorber_intervals_upper))
    vmr_ref           = read_field(ncid, 'vmr_ref', nlayers, nextabsorbers, ntemps)

    kmajor            = read_field(ncid, 'kmajor',  ngpts, nmixingfracs,  npress+1, ntemps)
    if(var_exists(ncid, 'rayl_lower')) then
      rayl_lower = read_field(ncid, 'rayl_lower',   ngpts, nmixingfracs,            ntemps)
      rayl_upper = read_field(ncid, 'rayl_upper',   ngpts, nmixingfracs,            ntemps)
    end if
    ! --------------------------------------------------
    !
    ! Initialize the gas optics class with data. The calls look slightly different depending
    !   on whether the radiation sources are internal to the atmosphere (longwave) or external (shortwave)
    ! gas_optics%load() returns a string; a non-empty string indicates an error.
    !
    if(var_exists(ncid, 'totplnk')) then
      !
      ! If there's a totplnk variable in the file it's a longwave (internal sources) type
      !
      totplnk     = read_field(ncid, 'totplnk', ninternalSourcetemps, nbnds)
      planck_frac = read_field(ncid, 'plank_fraction', ngpts, nmixingfracs, npress+1, ntemps)
      optimal_angle_fit = read_field(ncid, 'optimal_angle_fit', nfit_coeffs, nbnds)
      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor, &
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  totplnk, planck_frac,       &
                                  rayl_lower, rayl_upper, &
                                  optimal_angle_fit))
    else
      !
      ! Solar source doesn't have an dependencies yet
      !
      solar_quiet   = read_field(ncid, 'solar_source_quiet', ngpts)
      solar_facular = read_field(ncid, 'solar_source_facular', ngpts)
      solar_sunspot = read_field(ncid, 'solar_source_sunspot', ngpts)
      tsi_default   = read_field(ncid, 'tsi_default')
      mg_default    = read_field(ncid, 'mg_default')
      sb_default    = read_field(ncid, 'sb_default')
      call stop_on_err(kdist%load(available_gases, &
                                  gas_names,   &
                                  key_species, &
                                  band2gpt,    &
                                  band_lims,   &
                                  press_ref,   &
                                  press_ref_trop, &
                                  temp_ref,    &
                                  temp_ref_p, temp_ref_t,     &
                                  vmr_ref, kmajor,            &
                                  kminor_lower, kminor_upper, &
                                  gas_minor,identifier_minor,&
                                  minor_gases_lower, minor_gases_upper, &
                                  minor_limits_gpt_lower, &
                                  minor_limits_gpt_upper, &
                                  minor_scales_with_density_lower, &
                                  minor_scales_with_density_upper, &
                                  scaling_gas_lower, scaling_gas_upper, &
                                  scale_by_complement_lower, &
                                  scale_by_complement_upper, &
                                  kminor_start_lower, &
                                  kminor_start_upper, &
                                  solar_quiet, solar_facular, solar_sunspot, &
                                  tsi_default, mg_default, sb_default, &
                                  rayl_lower, rayl_upper))
    end if
    ! --------------------------------------------------
    ncid = nf90_close(ncid)
  end subroutine load_and_init_rrtmgp


!> @brief Read optical coefficients from NetCDF file for ecckd.
subroutine load_and_init_ecckd(ecckd, filename, available_gases)

  class(ty_gas_optics_ecckd), intent(inout) :: ecckd
  character(len=*), intent(in) :: filename
  class(ty_gas_concs), intent(in) :: available_gases

  integer :: band
  real(wp), dimension(:,:), allocatable :: band_lims_wvn
  integer, dimension(:,:), allocatable :: band2gpt
  character(len=1024) :: buffer
  character(len=32), dimension(16) :: composite_gas
  integer :: counter
  character(len=128) :: err
  logical :: found_gas
  character(len=32), dimension(16) :: gas
  integer, dimension(:), allocatable :: gpt2band
  integer :: i
  integer :: j
  integer, dimension(4) :: length
  logical :: lut
  integer :: n
  integer :: ncid
  integer :: num_tokens
  real, parameter :: total_solar_irradiance = 1361.d0
  logical :: uses_composite_gas
  character(len=64) :: varname

  if (nf90_open(trim(filename), NF90_NOWRITE, ncid) .ne. NF90_NOERR) then
    call stop_on_err("load_and_init_ecckd(): can't open file"//trim(fileName))
  endif

  length(1:1) = get_var_size(ncid, "pressure", 1)
  allocate(ecckd%log_pressure(length(1)))
  ecckd%log_pressure(:) = read_field(ncid, "pressure", length(1))
  ecckd%log_pressure = log(ecckd%log_pressure)

  length(1:2) = get_var_size(ncid, "temperature", 2)
  allocate(ecckd%temperature(length(1), length(2)))
  ecckd%temperature(:,:) = read_field(ncid, "temperature", length(1), length(2))

  length(1:1) = get_var_size(ncid, "wavenumber1_band", 1)
  allocate(band_lims_wvn(2, length(1)))
  band_lims_wvn(1,:) = read_field(ncid, "wavenumber1_band", length(1))
  band_lims_wvn(2,:) = read_field(ncid, "wavenumber2_band", length(1))
  allocate(band2gpt(2, length(1)))
  length(1:1) = get_var_size(ncid, "band_number", 1)
  allocate(gpt2band(length(1)))
  gpt2band(:) = read_field(ncid, "band_number", length(1))
  band2gpt(1, 1) = 1
  band2gpt(2, size(band2gpt, 2)) = length(1)
  band = 1
  do i = 2, length(1)
    if (gpt2band(i) .gt. band) then
      band2gpt(2, band) = i - 1
      band = band + 1
      band2gpt(1, band) = i
    endif
  enddo
  err = ecckd%init(band_lims_wvn, band2gpt)
  if (trim(err) .ne. "") then
    call stop_on_err(trim(err))
  endif
  deallocate(band_lims_wvn, gpt2band, band2gpt)

  length(1:2) = get_var_size(ncid, "gpoint_fraction", 2)
  allocate(ecckd%gpoint_fraction(length(1), length(2)))
  ecckd%gpoint_fraction(:,:) = read_field(ncid, "gpoint_fraction", length(1), length(2))

! if (variable_exists(dataset, "g_point")) then
!   allocate(ecckd%gpoint())
!   call read_data(dataset, "g_point", ecckd%gpoint)
!   allocate(ecckd%wavenumber_hr())
!   call read_data(dataset, "wavenumber_hr", ecckd%wavenumber_hr)
! endif

  ecckd%shortwave = var_exists(ncid, "solar_irradiance")
  if (ecckd%shortwave) then
    length(1:1) = get_var_size(ncid, "solar_irradiance", 1)
    allocate(ecckd%solar_irradiance(length(1)))
    ecckd%solar_irradiance(:) = read_field(ncid, "solar_irradiance", length(1))
    ecckd%total_solar_irradiance = sum(ecckd%solar_irradiance)
    length(1:1) = get_var_size(ncid, "rayleigh_molar_scattering_coeff", 1)
    allocate(ecckd%rayleigh_molar_scattering_coeff(length(1)))
    ecckd%rayleigh_molar_scattering_coeff(:) = read_field(ncid, "rayleigh_molar_scattering_coeff", &
                                                         length(1))
  else
    length(1:1) = get_var_size(ncid, "temperature_planck", 1)
    allocate(ecckd%temperature_planck(length(1)))
    ecckd%temperature_planck(:) = read_field(ncid, "temperature_planck", length(1))
    length(1:2) = get_var_size(ncid, "planck_function", 2)
    allocate(ecckd%planck_function(length(1), length(2)))
    ecckd%planck_function(:,:) = read_field(ncid, "planck_function", length(1), length(2))
  endif

  !Get the gases.
  buffer = ""
  call get_global_attribute(ncid, "constituent_id", buffer)
  gas(:) = ""
  call tokenize(buffer, gas, num_tokens)
  uses_composite_gas = .false.
  do i = 1, num_tokens
    if (trim(gas(i)) .eq. "composite") then
      uses_composite_gas = .true.
      call get_global_attribute(ncid, "composite_constituent_id", buffer)
      composite_gas(:) = ""
      call tokenize(buffer, composite_gas, ecckd%num_composite_gases)
      exit
    endif
  enddo
  counter = 0
  do i = 1, num_tokens
    if (trim(gas(i)) .ne. "composite") then
      counter = counter + 1
      ecckd%gas(counter) = trim(gas(i))
      call read_gas_input_data(ecckd%absorption(counter), ncid, gas(i))
      ecckd%absorption(counter)%composite_component = .false.
      ecckd%absorption(counter)%composite_only = .false.
      if (uses_composite_gas) then
        do j = 1, ecckd%num_composite_gases
          if (trim(gas(i)) .eq. trim(composite_gas(j))) then
            ecckd%absorption(counter)%composite_component = .true.
            exit
          endif
        enddo
      endif
    endif
  enddo
  if (uses_composite_gas) then
    do i = 1, ecckd%num_composite_gases
      found_gas = .false.
      do j = 1, num_tokens
        if (trim(composite_gas(i)) .eq. trim(gas(j))) then
          found_gas = .true.
          exit
        endif
      enddo
      if (.not. found_gas) then
        counter = counter + 1
        ecckd%gas(counter) = trim(composite_gas(i))
        call read_gas_input_data(ecckd%absorption(counter), ncid, "composite")
        ecckd%absorption(counter)%composite_component = .true.
        ecckd%absorption(counter)%composite_only = .true.
      endif
    enddo
  endif
  ecckd%num_gases = counter
  write(*,*) counter, ecckd%num_gases

  ncid = nf90_close(ncid)
end subroutine load_and_init_ecckd


subroutine read_gas_input_data(gas, ncid, gas_name)

  type(AbsorptionTable), intent(inout) :: gas
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: gas_name

  integer, dimension(4) :: length
  logical :: lut
  integer :: n
  character(len=132) :: varname

  lut = .false.
  varname = trim(gas_name)//"_mole_fraction"
  if (var_exists(ncid, varname)) then
    if (get_variable_num_dimensions(ncid, varname) .eq. 1) then
      lut = .true.
      gas%concentration_dependence_code = look_up_table
      varname = trim(gas_name)//"_mole_fraction"
      length(1:1) = get_var_size(ncid, varname, 1)
      allocate(gas%mole_fraction(length(1)))
      gas%mole_fraction(:) = read_field(ncid, varname, length(1))
      varname = trim(gas_name)//"_molar_absorption_coeff"
      length(1:4) = get_var_size(ncid, varname, 4)
      allocate(gas%coefficient(length(1), length(2), length(3), length(4)))
      gas%coefficient(:,:,:,:) = read_field(ncid, varname, length(1), &
                                            length(2), length(3), length(4))
    endif
  endif
  if (.not. lut) then
    varname = trim(gas_name)//"_conc_dependence_code"
    n = int(read_field(ncid, varname))
    if (n .eq. 0) then
      gas%concentration_dependence_code = none_
      varname = trim(gas_name)//"_mole_fraction"
      length(1:2) = get_var_size(ncid, varname, 2)
      allocate(gas%composite_vmr(length(1), length(2)))
      gas%composite_vmr(:,:) = read_field(ncid, varname, length(1), length(2))
    elseif (n .eq. 1) then
      gas%concentration_dependence_code = linear
    elseif (n .eq. 3) then
      gas%concentration_dependence_code = relative_linear
      varname = trim(gas_name)//"_reference_mole_fraction"
      gas%reference_mole_fraction = read_field(ncid, varname)
    else
      call stop_on_err("load_and_init_ecckd: bad concentration code for "//trim(gas_name))
    endif
    varname = trim(gas_name)//"_molar_absorption_coeff"
    n = get_variable_num_dimensions(ncid, varname)
    if (n .ne. 3) then
      call stop_on_err("load_and_init_ecckd: absorption coefficient not 3d for " &
                       //trim(gas_name))
    endif
    length(1:3) = get_var_size(ncid, varname, 3)
    allocate(gas%coefficient(length(1), length(2), length(3), 1))
    gas%coefficient(:,:,:,1) = read_field(ncid, varname, length(1), length(2), length(3))
  endif
end subroutine read_gas_input_data


!> @brief Read a global string attribute from a netCDF file.
subroutine get_global_attribute(ncid, attribute_name, attribute_value)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: attribute_name
  character(len=*), intent(inout) :: attribute_value

  integer :: err

  err = nf90_get_att(ncid, nf90_global, trim(attribute_name), attribute_value)
  if (err .ne. nf90_noerr) then
    call stop_on_err("get_global_attribute: error reading "//trim(attribute_name))
  endif
end subroutine get_global_attribute


function get_variable_num_dimensions(ncid, varname) result(n)

  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer :: n

  integer :: err
  integer :: varid

  err = nf90_inq_varid(ncid, trim(varname), varid)
  if (err .ne. nf90_noerr) then
    call stop_on_err("get_variable_num_dimensions: error finding variable "//trim(varname))
  endif
  err = nf90_inquire_variable(ncid, varid, ndims=n)
  if (err .ne. nf90_noerr) then
    call stop_on_err("get_variable_num_dimensions: error finding variable "//trim(varname))
  endif
end function get_variable_num_dimensions


!> @brief Break a buffer in to tokens.
subroutine tokenize(buffer, tokens, num_tokens)

  character(len=*), intent(in) :: buffer !< String buffer.
  character(len=*), dimension(:), intent(inout) :: tokens !< Array of string tokens.
  integer, intent(out) :: num_tokens

  integer :: counter
  logical :: found_token
  integer :: i
  integer :: n
  integer :: start

  n = len_trim(buffer)
  found_token = .false.
  start = 1
  counter = 0
  do i = 1, n
    if (found_token) then
      if (i .eq. n .and. buffer(i:i) .ne. " ") then
        ! Last character is not a space.
        counter = counter + 1
        if (counter .gt. size(tokens)) then
          call stop_on_err("tokenize: more tokens than can fit in array.")
        endif
        if (i - start + 1 .gt. len(tokens(counter))) then
          call stop_on_err("tokenize: token is cut off.")
        endif
        tokens(counter) = buffer(start:i)
        exit
      elseif (buffer(i:i) .eq. " ") then
        ! Reached the end of a token.
        counter = counter + 1
        if (counter .gt. size(tokens)) then
          call stop_on_err("tokenize: more tokens than can fit in array.")
        endif
        if (i - start .gt. len(tokens(counter))) then
          call stop_on_err("tokenize: token is cut off.")
        endif
        tokens(counter) = buffer(start:i-1)
        found_token = .false.
      endif
    else
      if (buffer(i:i) .ne. " ") then
        found_token = .true.
        start = i
      endif
    endif
  enddo
  num_tokens = counter
end subroutine tokenize


end module

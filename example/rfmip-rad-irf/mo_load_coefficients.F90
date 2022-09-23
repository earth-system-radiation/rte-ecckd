module mo_load_coefficients
use mo_rte_kind, only: wp
use mo_gas_concentrations, only: ty_gas_concs
use gas_optics_ecckd, only: AbsorptionTable, none_, linear, look_up_table, &
                            relative_linear, ty_gas_optics_ecckd
use mo_simple_netcdf, only: read_field, read_field_int, var_exists, get_var_size, stop_on_err
use netcdf
implicit none
private


public :: load_and_init


contains


!> @brief Read optical coefficients from NetCDF file for ecckd.
subroutine load_and_init(ecckd, filename, available_gases)

  class(ty_gas_optics_ecckd), intent(inout) :: ecckd
  character(len=*), intent(in) :: filename !< Path to input file.
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
  integer :: ncid
  integer :: num_tokens
  logical :: uses_composite_gas

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
  gpt2band(:) = read_field_int(ncid, "band_number", length(1))
  gpt2band(:) = gpt2band(:) + 1
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
end subroutine load_and_init


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
!     allocate(gas%composite_vmr(length(1), length(2)))
!     gas%composite_vmr(:,:) = read_field(ncid, varname, length(1), length(2))
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


!> @brief Get the number of dimensions for a variable in netCDF dataset.
function get_variable_num_dimensions(ncid, varname) result(n)

  integer, intent(in) :: ncid !< NetCDF dataset id.
  character(len=*), intent(in) :: varname !< Variable name.
  integer :: n !< Number of dimensions.

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
  integer, intent(out) :: num_tokens !< Number of tokens found.

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


end module mo_load_coefficients

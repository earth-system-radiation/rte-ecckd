module gas_optics_ecckd
use, intrinsic :: iso_fortran_env, only: error_unit
use mo_gas_concentrations, only: ty_gas_concs
use mo_gas_optics, only: ty_gas_optics
use mo_optical_props, only: ty_optical_props_arry, ty_optical_props_2str
use mo_rte_kind, only: wp
use mo_source_functions, only: ty_source_func_lw
implicit none
private


!> @brief Helper class used to calculate optical depth for each gas.
type, public :: AbsorptionTable
  real(kind=wp), dimension(:,:,:,:), allocatable :: coefficient !< Absorption coefficient [m2 mol-1] (gpoint, pressure, temperature).
  logical :: composite_only !< Only part of the "composite" gas?
  integer :: concentration_dependence_code !< How is absoprtion coefficient calculated?
  real(kind=wp), dimension(:), allocatable :: mole_fraction !< Mole fraction look up table.
  real(kind=wp) :: reference_mole_fraction !< Reference mole fraction amount.
end type AbsorptionTable


!> @brief Class implementing the ecckd correlated-k distribution gas optics model.
type, extends(ty_gas_optics), public :: ty_gas_optics_ecckd
  type(AbsorptionTable), dimension(16) :: absorption !< Individual gas data.
  character(len=32), dimension(16) :: gas !< Gas names.
  real(kind=wp), dimension(:,:), allocatable :: gpoint_fraction
  real(kind=wp), dimension(:), allocatable :: log_pressure !< Natural log of pressure [Pa] (pressure).
  integer :: num_composite_gases !< Number of gases that make up the "composite" gas.
  integer :: num_gases !< Total number of gases.
  real(kind=wp), dimension(:,:), allocatable :: planck_function !< Planck function flux into a horizontal plane [W m-2] (gpoint, planck temperature).
  real(kind=wp), dimension(:), allocatable :: rayleigh_molar_scattering_coeff !< Rayleigh molar scattering coefficient [m2 mol-1] (gpoint).
  logical :: shortwave !< Longwave or shortwave.
  real(kind=wp), dimension(:), allocatable :: solar_irradiance !< Solar irradiance [W m-2] (gpoint).
  real(kind=wp), dimension(:,:), allocatable :: temperature !< Temperature [K] (pressure, temperature).
  real(kind=wp), dimension(:), allocatable :: temperature_planck !< Planck temperature [K] (planck temperature).
  real(kind=wp) :: total_solar_irradiance !< Total solar irradiance [W m-2].
  contains
  procedure, public :: source_is_internal
  procedure, public :: source_is_external
  procedure, public :: get_ngas
  procedure, public :: get_gases
  procedure, public :: get_press_min
  procedure, public :: get_press_max
  procedure, public :: get_temp_min
  procedure, public :: get_temp_max
  procedure, public :: gas_optics_int
  procedure, public :: gas_optics_ext
end type ty_gas_optics_ecckd


real(kind=wp), parameter :: gravity = 9.80665 !< Acceleration due to gravity [m s-2].
real(kind=wp), parameter :: dry_air_molar_mass = 28.970 !< Dry air molar mass [g mol-1].
real(kind=wp), parameter :: pi = 3.14159265359
integer, parameter, public :: none_ = 0
integer, parameter, public :: linear = 1
integer, parameter, public :: look_up_table = 2
integer, parameter, public :: relative_linear = 3


contains


!> @brief Calculate the optical depth for a gas at each g-point.
subroutine calculate_optical_depth(this, gas, level_pressure, layer_temperature, &
                                   layer_vmr, optical_depth, logarithmic_interpolation)

  type(ty_gas_optics_ecckd), intent(in) :: this
  integer, intent(in) :: gas !< Index of gas to calculate.
  real(kind=wp), dimension(:,:), intent(in) :: level_pressure !< Pressure [Pa] (column, level).
  real(kind=wp), dimension(:,:), intent(in) :: layer_temperature !< Temperature [K] (column, layer).
  real(kind=wp), dimension(:,:), intent(in) :: layer_vmr !< Volume-mixing ratio [mol mol-1] (column, layer).
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: optical_depth !< Optical depth (column, layer, gpoint).
  logical, intent(in) :: logarithmic_interpolation !< Use logarithmic interpolation.

  real(kind=wp) :: d_log_p
  real(kind=wp) :: d_log_vmr
  real(kind=wp) :: dt
  real(kind=wp) :: global_weight
  integer :: i
  integer :: ip0
  integer :: it0
  integer :: iv0
  integer :: j
  integer :: k
  real(kind=wp) :: log_p_0
  real(kind=wp) :: log_pressure
  real(kind=wp) :: log_vmr
  integer :: num_columns
  integer :: num_gpoints
  integer :: num_layers
  real(kind=wp) :: pressure_index
  real(kind=wp) :: pressure_weight0
  real(kind=wp) :: pressure_weight1
  real(kind=wp) :: simple_weight
  real(kind=wp) :: temperature_index
  real(kind=wp) :: temperature_weight0
  real(kind=wp) :: temperature_weight1
  real(kind=wp) :: t0
  real(kind=wp) :: vmr_index
  real(kind=wp) :: vmr_weight0
  real(kind=wp) :: vmr_weight1
  real(kind=wp) :: weight

  log_p_0 = this%log_pressure(1)
  d_log_p = this%log_pressure(2) - this%log_pressure(1)
  dt = this%temperature(1,2) - this%temperature(1,1)
  global_weight = 1./(gravity*0.001*dry_air_molar_mass)

  num_columns = size(layer_temperature, 1)
  num_gpoints = size(this%gpoint_fraction, 2)
  num_layers = size(layer_temperature, 2)
  if (allocated(optical_depth)) then
    deallocate(optical_depth)
  endif
  allocate(optical_depth(num_columns, num_layers, num_gpoints))

  do j = 1, num_layers
    do i = 1, num_columns
      ! Pressure interpolation points.
      log_pressure = log(0.5*(level_pressure(i,j+1) + level_pressure(i,j)))
      pressure_index = (log_pressure - log_p_0)/d_log_p
      pressure_index = 1. + &
                         max(0._wp, &
                             min(pressure_index, &
                                 real(size(this%log_pressure) - 1.0001_wp, kind=wp)))
      ip0 = int(pressure_index)
      pressure_weight1 = pressure_index - ip0
      pressure_weight0 = 1. - pressure_weight1

      ! Temperature interpolation points.
      t0 = pressure_weight0*this%temperature(ip0,1) + &
           pressure_weight1*this%temperature(ip0+1,1)
      temperature_index = (layer_temperature(i,j) - t0)/dt
      temperature_index = 1. + &
                            max(0._wp, &
                              min(temperature_index, &
                                  real(size(this%temperature, 2) - 1.0001_wp, kind=wp)))
      it0 = int(temperature_index)
      temperature_weight1 = temperature_index - it0
      temperature_weight0 = 1. - temperature_weight1

      ! Weighting.
      simple_weight = global_weight*(level_pressure(i,j+1) - level_pressure(i,j))
      weight = 0.
      if (this%absorption(gas)%concentration_dependence_code .eq. relative_linear) then
        weight = simple_weight*(layer_vmr(i,j) - this%absorption(gas)%reference_mole_fraction)
      else
        weight = simple_weight*layer_vmr(i,j)
      endif

      if (this%absorption(gas)%concentration_dependence_code .eq. look_up_table) then
        ! Volume mixing ratio interpolation points.
        log_vmr = log(max(layer_vmr(i,j), this%absorption(gas)%mole_fraction(1)))
        d_log_vmr = log(this%absorption(gas)%mole_fraction(2)/ &
                        this%absorption(gas)%mole_fraction(1))
        vmr_index = (log_vmr - log(this%absorption(gas)%mole_fraction(1)))/d_log_vmr
        vmr_index = 1. + &
                      max(0._wp, &
                        min(vmr_index, &
                            real(size(this%absorption(gas)%mole_fraction) - 1.001_wp, kind=wp)))
        iv0 = int(vmr_index)
        vmr_weight1 = vmr_index - iv0
        vmr_weight0 = 1. - vmr_weight1

        !Tri-linear absorption interpolation.
        if (.not. logarithmic_interpolation) then
          optical_depth(i,j,:) = weight*(vmr_weight0*(temperature_weight0*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0,iv0) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0,iv0)) + &
            temperature_weight1*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0+1,iv0) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0+1,iv0))) + &
            vmr_weight1*(temperature_weight0*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0,iv0+1) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0,iv0+1)) + &
            temperature_weight1*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0+1,iv0+1) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0+1,iv0+1))))
        else
          optical_depth(i,j,:) = weight* &
            exp(vmr_weight0*(temperature_weight0*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0,iv0)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0,iv0))) + &
            temperature_weight1*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0+1,iv0)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0+1,iv0)))) + &
            vmr_weight1*(temperature_weight0*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0,iv0+1)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0,iv0+1))) + &
            temperature_weight1*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0+1,iv0+1)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0+1,iv0+1)))))
        endif
      elseif (this%absorption(gas)%concentration_dependence_code .eq. linear .or. &
              this%absorption(gas)%concentration_dependence_code .eq. relative_linear) then
        if (.not. logarithmic_interpolation) then
          !Bi-linear interpolation.
          optical_depth(i,j,:) = weight*(temperature_weight0*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0,1) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0,1)) + &
            temperature_weight1*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0+1,1) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0+1,1)))
        else
          optical_depth(i,j,:) = weight* &
            exp(temperature_weight0*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0,1)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0,1))) + &
            temperature_weight1*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0+1,1)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0+1,1))))
        endif
      else
        if (.not. logarithmic_interpolation) then
          !Bi-linear interpolation.
          optical_depth(i,j,:) = simple_weight*(temperature_weight0*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0,1) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0,1)) + &
            temperature_weight1*( &
            pressure_weight0*this%absorption(gas)%coefficient(:,ip0,it0+1,1) + &
            pressure_weight1*this%absorption(gas)%coefficient(:,ip0+1,it0+1,1)))
        else
          optical_depth(i,j,:) = simple_weight* &
            exp(temperature_weight0*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0,1)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0,1))) + &
            temperature_weight1*( &
            pressure_weight0*log(this%absorption(gas)%coefficient(:,ip0,it0+1,1)) + &
            pressure_weight1*log(this%absorption(gas)%coefficient(:,ip0+1,it0+1,1))))
        endif
      endif

      !Remove negative optical depths.
      do k = 1, size(optical_depth, 3)
        if (optical_depth(i,j,k) .lt. 0.) then
          optical_depth(i,j,k) = 0.
        endif
      enddo
    enddo
  enddo
end subroutine calculate_optical_depth


!> @brief Calculate the planck function.
subroutine calculate_planck_function(this, level_temperature, planck)

  type(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp), dimension(:,:), intent(in) :: level_temperature !< Temperature [K] (column, level).
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: planck !< Intensity [W m-2 sr-1] (column, level, gpoint).

  real(kind=wp) :: dt
  integer :: i
  integer :: it0
  integer :: j
  integer :: num_columns
  integer :: num_gpoints
  integer :: num_levels
  real(kind=wp) :: temperature_index
  real(kind=wp) :: temperature_weight0
  real(kind=wp) :: temperature_weight1
  real(kind=wp) :: t0

  num_columns = size(level_temperature, 1)
  num_levels = size(level_temperature, 2)
  num_gpoints = size(this%gpoint_fraction, 2)
  if (allocated(planck)) then
    deallocate(planck)
  endif
  allocate(planck(num_columns, num_levels, num_gpoints))

  dt = this%temperature_planck(2) - this%temperature_planck(1)
  t0 = this%temperature_planck(1)
  do j = 1, num_levels
    do i = 1, num_columns
      temperature_index = (level_temperature(i,j) - t0)/dt
      if (temperature_index .ge. 0) then
        temperature_index = 1. + temperature_index
        it0 = min(int(temperature_index), size(this%temperature_planck) - 1)
        temperature_weight1 = temperature_index - it0
        temperature_weight0 = 1. - temperature_weight1
        planck(i,j,:) = temperature_weight0*this%planck_function(:,it0) + &
                        temperature_weight1*this%planck_function(:,it0+1)
      else
        planck(i,j,:) = (level_temperature(i,j)/t0)*this%planck_function(:,1)
      endif
    enddo
  enddo
  planck(:,:,:) = planck(:,:,:)/pi
end subroutine calculate_planck_function


!> @brief Calculate the rayleigh scattering optical depth.
subroutine calculate_rayleigh_optical_depth(this, level_pressure, optical_depth)

  type(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp), dimension(:,:), intent(in) :: level_pressure !< Pressure [Pa] (column, level).
  real(kind=wp), dimension(:,:,:), allocatable, intent(inout) :: optical_depth !< Optical depth (column, layer, gpoint).

  integer :: i
  real(kind=wp), dimension(:,:), allocatable :: moles_per_layer
  integer :: num_columns
  integer :: num_gpoints
  integer :: num_layers

  num_columns = size(level_pressure, 1)
  num_layers = size(level_pressure, 2) - 1
  num_gpoints = size(this%gpoint_fraction, 2)
  if (allocated(optical_depth)) then
    deallocate(optical_depth)
  endif
  allocate(optical_depth(num_columns, num_layers, num_gpoints))
  allocate(moles_per_layer(num_columns, num_layers))
  moles_per_layer = (level_pressure(:,2:num_layers+1) - level_pressure(:,1:num_layers))* &
                    (1./(gravity*0.001*dry_air_molar_mass))
  do i = 1, num_gpoints
    optical_depth(:,:,i) = moles_per_layer(:,:)*this%rayleigh_molar_scattering_coeff(i)
  enddo
  deallocate(moles_per_layer)
end subroutine calculate_rayleigh_optical_depth


!> @brief Calculate the optical depth for all gases.
function gas_optical_depth(this, plev, tlay, gas_desc, optical_props) &
  result(error_msg)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp), dimension(:,:), intent(in) :: plev !< Level pressures [Pa]; (ncol, nlay+1)
  real(kind=wp), dimension(:,:), intent(in) :: tlay !< Layer temperatures [K]; (ncol, nlay)
  type(ty_gas_concs), intent(in) :: gas_desc !< Gas volume mixing ratios
  class(ty_optical_props_arry), intent(inout) :: optical_props !< Optical properties
  character(len=128) :: error_msg !< Error string (empty if succssful)

  logical :: first_calc
  real(kind=wp), dimension(:,:), allocatable :: layer_vmr
  integer :: i, j, n
  character(len=32), dimension(:), allocatable :: names
  real(kind=wp), dimension(:,:,:), allocatable :: optical_depth

  !Calculate the gas optical depth.
  n = gas_desc%get_num_gases()
  allocate(names(n))
  names = gas_desc%get_gas_names()
  allocate(layer_vmr(size(tlay, 1), size(tlay, 2)))
  allocate(optical_depth(size(tlay, 1), size(tlay, 2), size(this%gpoint_fraction, 2)))
  optical_depth(:,:,:) = 0.
  optical_props%tau(:,:,:) = 0.
  first_calc = .true.
  do j = 1, n
    do i = 1, this%num_gases
      if (trim(this%gas(i)) .eq. trim(names(j))) then
        error_msg = gas_desc%get_vmr(names(j), layer_vmr)
        if (trim(error_msg) .ne. "") then
          return
        endif
        exit
      endif
    enddo
    if (i .gt. this%num_gases) then
      !Did not find the gas
!     write(error_unit, *) "Warning: Failed to find "//trim(names(j))//" in the ecckd model"
      cycle
!     error_msg = "Failed to find "//trim(names(j))//" in the ecckd model"
!     return
    endif
    if (this%absorption(i)%composite_only .and. .not. first_calc) then
      cycle
    endif
    call calculate_optical_depth(this, i, plev, tlay, &
                                 layer_vmr, optical_depth, .false.)
    optical_props%tau(:,:,:) = optical_props%tau(:,:,:) + optical_depth(:,:,:)
    if (this%absorption(i)%composite_only) then
      first_calc = .false.
    endif
  enddo
  deallocate(names, layer_vmr, optical_depth)
end function gas_optical_depth


!> @brief Compute gas optical depth and Planck source functions,
!!        given temperature, pressure, and composition.
function gas_optics_int(this, play, plev, tlay, tsfc, gas_desc, &
                        optical_props, sources, col_dry, tlev) &
  result(error_msg)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp), dimension(:,:), intent(in) :: play !< Layer pressures [Pa]; (ncol, nlay)
  real(kind=wp), dimension(:,:), intent(in) :: plev !< Level pressures [Pa]; (ncol, nlay+1)
  real(kind=wp), dimension(:,:), intent(in) :: tlay !< Layer temperatures [K]; (ncol, nlay)
  real(kind=wp), dimension(:), intent(in) :: tsfc !< Surface skin temperatures [K]; (ncol)
  type(ty_gas_concs), intent(in) :: gas_desc !< Gas volume mixing ratios
  class(ty_optical_props_arry), intent(inout) :: optical_props !< Optical properties
  class(ty_source_func_lw), intent(inout) :: sources !< Planck sources
  character(len=128) :: error_msg !< Error string (empty if succssful)
  real(kind=wp), dimension(:,:), intent(in), target, optional :: col_dry !< Column dry amount [cm-2]; (col, nlay)
  real(kind=wp), dimension(:,:), intent(in), target, optional :: tlev !< Level temperatures [K]; (ncol, nlay+1)

  real(kind=wp), dimension(:,:,:), allocatable :: buffer
  real(kind=wp), dimension(:,:), allocatable :: surface_temperature

  !Calculate optical depth.
  error_msg = gas_optical_depth(this, plev, tlay, gas_desc, optical_props)
  if (len_trim(error_msg) .gt. 0) then
    return
  endif

  !Calculate source functions.
  call calculate_planck_function(this, tlay, sources%lay_source)
  allocate(surface_temperature(size(tsfc), 1))
  surface_temperature(:,1) = tsfc(:)
  allocate(buffer(size(sources%sfc_source, 1), 1, size(sources%sfc_source, 2)))
  call calculate_planck_function(this, surface_temperature, buffer)
  sources%sfc_source(:,:) = buffer(:,1,:)
  deallocate(buffer, surface_temperature)
  if (.not. present(tlev)) then
    error_msg = "tlev is required for ecckd"
    return
  endif

  allocate(buffer(size(sources%lev_source_inc, 1), &
                  size(sources%lev_source_inc, 2) + 1, &
                  size(sources%lev_source_inc, 3)))
  call calculate_planck_function(this, tlev, buffer)
  sources%lev_source_inc(:,:,:) = buffer(:,2:,:)
  sources%lev_source_dec(:,:,:) = buffer(:,:size(sources%lev_source_dec, 2),:)
  deallocate(buffer)
end function gas_optics_int


!> @brief Compute gas optical depth given temperature, pressure, and composition.
!!        Top-of-atmosphere stellar insolation is also reported.
function gas_optics_ext(this, play, plev, tlay, gas_desc, optical_props, toa_src, col_dry) &
  result(error_msg)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp), dimension(:,:), intent(in) :: play !< Layer pressures [Pa]; (ncol, nlay)
  real(kind=wp), dimension(:,:), intent(in) :: plev !< Level pressures [Pa]; (ncol, nlay+1)
  real(kind=wp), dimension(:,:), intent(in) :: tlay !< Layer temperatures [K]; (ncol, nlay)
  type(ty_gas_concs), intent(in) :: gas_desc !< Gas volume mixing ratios
  class(ty_optical_props_arry), intent(inout) :: optical_props
  real(kind=wp), dimension(:,:), intent(out) :: toa_src !< Incoming solar irradiance (ncol, ngpt)
  real(kind=wp), dimension(:,:), intent(in), target, optional :: col_dry !< Column dry amount [cm-2]; (col, nlay)
  character(len=128) :: error_msg !< String error message (empty if successful)

  integer :: i
  integer :: j
  real(kind=wp), dimension(:,:,:), allocatable :: optical_depth

  !Calculate optical depth.
  error_msg = gas_optical_depth(this, plev, tlay, gas_desc, optical_props)
  if (len_trim(error_msg) .gt. 0) then
    return
  endif

  !Calculate the rayleigh optical depth and single scatter albedo.
  call calculate_rayleigh_optical_depth(this, plev, optical_depth)
  optical_props%tau(:,:,:) = optical_props%tau(:,:,:) + optical_depth(:,:,:)
  select type(optical_props)
    type is (ty_optical_props_2str)
      optical_props%ssa(:,:,:) = optical_depth(:,:,:)/optical_props%tau(:,:,:)
      optical_props%g(:,:,:) = 0
    class default
      error_msg = "shortwave must use ty_optical_props_2str"
      return
  end select
  deallocate(optical_depth)

  !Calculate TOA source function.
  do j = 1, size(this%gpoint_fraction, 2)
    do i = 1, size(tlay, 1)
      toa_src(i,j) = this%solar_irradiance(j)
    enddo
  enddo
end function gas_optics_ext


!> @brief Return the number of gases registered in the spectral configuration.
pure function get_ngas(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  integer :: get_ngas

  get_ngas = this%num_gases
end function get_ngas


!> @brief Return true if initialized for internal sources/longwave, false otherwise.
pure function source_is_internal(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  logical :: source_is_internal

  source_is_internal = allocated(this%temperature_planck)
end function source_is_internal


!> @brief Return true if initialized for external sources/shortwave, false otherwise.
pure function source_is_external(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  logical :: source_is_external

  source_is_external = allocated(this%solar_irradiance)
end function source_is_external


!> @brief Return the names of the gases known to the k-distributions.
pure function get_gases(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  character(len=32), dimension(this%num_gases) :: get_gases

  get_gases = this%gas(1:this%num_gases)
end function get_gases


!> @brief Return the minimum pressure on the interpolation grids.
pure function get_press_min(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp) :: get_press_min

  get_press_min = exp(this%log_pressure(1))
end function get_press_min


!> @brief Return the maximum pressure on the interpolation grids.
pure function get_press_max(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp) :: get_press_max

  get_press_max = exp(this%log_pressure(size(this%log_pressure)))
end function get_press_max


!> @brief Return the minimum temparature on the interpolation grids.
pure function get_temp_min(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp) :: get_temp_min

  get_temp_min = minval(this%temperature)
end function get_temp_min


!> @brief Return the maximum temparature on the interpolation grids.
pure function get_temp_max(this)

  class(ty_gas_optics_ecckd), intent(in) :: this
  real(kind=wp) :: get_temp_max

  get_temp_max = maxval(this%temperature)
end function get_temp_max


end module gas_optics_ecckd

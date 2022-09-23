module utils
use mo_simple_netcdf, only: stop_on_err
implicit none
private


public :: determine_gas_names
public :: parse_args


contains


!> @brief Print usage instructions.
subroutine usage()
  use, intrinsic :: iso_fortran_env, only : error_unit

  character(len=32) :: command

  call get_command_argument(0, command)
  write(error_unit, *) "Usage: "//trim(command)//" rfmip_file ecckd_file"
end subroutine usage


!> @brief Print help message.
subroutine help()
  use, intrinsic :: iso_fortran_env, only : error_unit

  call usage()
  write(error_unit, *) ""
  write(error_unit, *) "Args:"
  write(error_unit, *) "rfmip_file - RFMIP input file."
  write(error_unit, *) "ecckd_file - ecckd input file."
  write(error_unit, *) "-f [1,2] - Forcing index."
  write(error_unit, *) "-h|--help - Prints this help message."
  write(error_unit, *) "-p [1,2] - Physics index."
end subroutine help


!> @brief Set gas names to look for in the input files.
subroutine determine_gas_names(forcing_index, names_in_kdist, names_in_rfmip)

  integer, intent(in) :: forcing_index !< Forcing index [1,2].
  character(len=32), dimension(:), intent(inout) :: names_in_kdist !< Gas names in the code.
  character(len=32), dimension(:), intent(inout) :: names_in_rfmip !< Gas names in the RFMIP file.

  names_in_kdist = (/"co2  ", &
                     "ch4  ", &
                     "n2o  ", &
                     "o2   ", &
                     "cfc11", &
                     "cfc12"/)
  if (forcing_index .eq. 1) then
    names_in_rfmip = (/"carbon_dioxide", &
                       "methane       ", &
                       "nitrous_oxide ", &
                       "oxygen        ", &
                       "cfc11         ", &
                       "cfc12         "/)
  elseif (forcing_index .eq. 2) then
    names_in_rfmip = (/"carbon_dioxide", &
                       "methane       ", &
                       "nitrous_oxide ", &
                       "oxygen        ", &
                       "cfc11eq       ", &
                       "cfc12         "/)
  else
    call stop_on_err("forcing index must equal 1 or 2.")
  endif
end subroutine determine_gas_names


!> @brief Parses command line arguments.
subroutine parse_args(rfmip_path, ecckd_path, forcing_index, physics_index)

  character(len=*), intent(inout) :: rfmip_path !< Path to the RFMIP input file.
  character(len=*), intent(inout) :: ecckd_path !< Path to the ecckd input file.
  integer, intent(inout) :: forcing_index !< Forcing index.
  integer, intent(inout) :: physics_index !< Physics index.

  character(len=128) :: buffer
  integer :: i, j
  integer :: nargs

  forcing_index = 1
  physics_index = 1
  nargs = command_argument_count()
  if (nargs .lt. 2) then
    call usage()
    stop 1
  endif
  i = 1
  j = 0
  do while (i .le. nargs)
    call get_command_argument(i, buffer)
    if (trim(buffer) .eq. "-h" .or. trim(buffer) .eq. "--help") then
      call help()
      stop 0
    elseif (trim(buffer) .eq. "-f") then
      i = i + 1
      if (i .gt. nargs) then
        call usage()
        stop 1
      endif
      call get_command_argument(i, buffer)
      read(buffer, "(i1)") forcing_index
      if (forcing_index .lt. 1 .or. forcing_index .gt. 2) then
        call stop_on_err("forcing index must be either 1 or 2.")
      endif
    elseif (trim(buffer) .eq. "-p") then
      i = i + 1
      if (i .gt. nargs) then
        call usage()
        stop 1
      endif
      call get_command_argument(i, buffer)
      read(buffer, "(i1)") physics_index
      if (physics_index .lt. 1 .or. physics_index .gt. 2) then
        call stop_on_err("physics index must be either 1 or 2.")
      endif
    else
      if (j .eq. 0) then
        rfmip_path = trim(buffer)
      elseif (j .eq. 1) then
        ecckd_path = trim(buffer)
      else
        call usage()
        stop 1
      endif
      j = j + 1
    endif
    i = i + 1
  enddo
end subroutine parse_args


end module utils

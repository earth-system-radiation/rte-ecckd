module simple_netcdf
use mo_rte_kind, only: wp
use netcdf
implicit none
private


interface read_field
  module procedure read_scalar
  module procedure read_1d_field
  module procedure read_2d_field
  module procedure read_3d_field
  module procedure read_4d_field
end interface


interface read_field_int
  module procedure read_1d_field_int
end interface


interface write_field
  module procedure write_1d_int_field
  module procedure write_2d_int_field
  module procedure write_1d_field
  module procedure write_2d_field
  module procedure write_3d_field
  module procedure write_4d_field
end interface


public :: get_dim_size
public :: get_var_size
public :: read_field
public :: read_field_int
public :: stop_on_err
public :: var_exists
public :: write_field


contains


function read_scalar(ncid, varname)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real(kind=wp) :: read_scalar

  integer :: varid

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't find variable "//trim(varname))
  endif
  if (nf90_get_var(ncid, varid, read_scalar) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't read variable "//trim(varname))
  endif
end function read_scalar


function read_1d_field(ncid, varname, nx)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nx
  real(kind=wp), dimension(nx) :: read_1d_field

  integer :: varid

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't find variable "//trim(varname))
  endif
  if (nf90_get_var(ncid, varid, read_1d_field) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't read variable "//trim(varname))
  endif
end function read_1d_field


function read_1d_field_int(ncid, varname, nx)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nx
  integer, dimension(nx) :: read_1d_field_int

  integer :: varid

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't find variable "//trim(varname))
  endif
  if (nf90_get_var(ncid, varid, read_1d_field_int) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't read variable "//trim(varname))
  endif
end function read_1d_field_int


function read_2d_field(ncid, varname, nx, ny)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nx, ny
  real(kind=wp), dimension(nx, ny) :: read_2d_field

  integer :: varid

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't find variable "//trim(varname))
  endif
  if (nf90_get_var(ncid, varid, read_2d_field) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't read variable "//trim(varname))
  endif
end function read_2d_field


function read_3d_field(ncid, varname, nx, ny, nz)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nx, ny, nz
  real(kind=wp), dimension(nx, ny, nz) :: read_3d_field

  integer :: varid

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't find variable "//trim(varname))
  endif
  if (nf90_get_var(ncid, varid, read_3d_field) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't read variable "//trim(varname))
  endif
end function read_3d_field


function read_4d_field(ncid, varname, nw, nx, ny, nz)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: nw, nx, ny, nz
  real(kind=wp), dimension(nw, nx, ny, nz) :: read_4d_field

  integer :: varid

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't find variable "//trim(varname))
  endif
  if (nf90_get_var(ncid, varid, read_4d_field) .ne. nf90_noerr) then
    call stop_on_err("read_field: can't read variable "//trim(varname))
  endif
end function read_4d_field


function write_1d_int_field(ncid, varname, var) result(err_msg)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, dimension(:), intent(in) :: var
  character(len=128) :: err_msg

  integer :: varid

  err_msg = ""
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    err_msg = "write_field: can't find variable "//trim(varname)
    return
  endif
  if (nf90_put_var(ncid, varid, var) .ne. nf90_noerr) then
    err_msg = "write_field: can't write variable "//trim(varname)
  endif
end function write_1d_int_field


function write_2d_int_field(ncid, varname, var) result(err_msg)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, dimension(:,:), intent(in) :: var
  character(len=128) :: err_msg

  integer :: varid

  err_msg = ""
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    err_msg = "write_field: can't find variable "//trim(varname)
    return
  endif
  if (nf90_put_var(ncid, varid, var) .ne. nf90_noerr) then
    err_msg = "write_field: can't write variable "//trim(varname)
  endif
  end function write_2d_int_field


function write_1d_field(ncid, varname, var) result(err_msg)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real(kind=wp), dimension(:), intent(in) :: var
  character(len=128) :: err_msg

  integer :: varid

  err_msg = ""
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    err_msg = "write_field: can't find variable "//trim(varname)
    return
  endif
  if (nf90_put_var(ncid, varid, var) .ne. nf90_noerr) then
    err_msg = "write_field: can't write variable "//trim(varname)
  endif
end function write_1d_field


function write_2d_field(ncid, varname, var) result(err_msg)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real(kind=wp), dimension(:,:), intent(in) :: var
  character(len=128) :: err_msg

  integer :: varid
  integer :: stat

  err_msg = ""
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    err_msg = "write_field: can't find variable "//trim(varname)
    return
  endif
  stat = nf90_put_var(ncid, varid, var)
  if (stat .ne. nf90_noerr) then
    err_msg = "write_field: can't write variable "//trim(varname)// &
              " netcdf err: "//nf90_strerror(stat)
  endif
end function write_2d_field


function write_3d_field(ncid, varname, var) result(err_msg)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real(kind=wp), dimension(:,:,:), intent(in) :: var
  character(len=128) :: err_msg

  integer :: varid
  integer :: stat

  err_msg = ""
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    err_msg = "write_field: can't find variable "//trim(varname)
    return
  endif
  stat = nf90_put_var(ncid, varid, var)
  if (stat .ne. nf90_noerr) then
    err_msg = "write_field: can't write variable "//trim(varname)// &
              " netcdf err: "//nf90_strerror(stat)
  endif
end function write_3d_field


function write_4d_field(ncid, varname, var) result(err_msg)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real(kind=wp), dimension(:,:,:,:), intent(in) :: var
  character(len=128) :: err_msg

  integer :: varid
  integer :: stat

  err_msg = ""
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    err_msg = "write_field: can't find variable "//trim(varname)
    return
  endif
  stat = nf90_put_var(ncid, varid, var)
  if (stat .ne. nf90_noerr) then
    err_msg = "write_field: can't write variable "//trim(varname)// &
              " netcdf err: "//nf90_strerror(stat)
  endif
end function write_4d_field


!> @brief does this variable exist (have a valid var_id) in the open netcdf file?
function var_exists(ncid, varname)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  logical :: var_exists

  integer :: varid

  var_exists = nf90_inq_varid(ncid, trim(varname), varid) .eq. nf90_noerr
end function var_exists


!> @brief get the length of a dimension from an open netcdf file.
function get_dim_size(ncid, dimname)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: dimname
  integer :: get_dim_size

  integer :: dimid

  if (nf90_inq_dimid(ncid, trim(dimname), dimid) .eq. nf90_noerr) then
    if (nf90_inquire_dimension(ncid, dimid, len=get_dim_size) .ne. nf90_noerr) then
      get_dim_size = 0
    endif
  else
    get_dim_size = 0
  endif
end function get_dim_size


!> @brief returns the extents of a netcdf variable on disk.
function get_var_size(ncid, varname, n)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  integer, intent(in) :: n
  integer, dimension(n) :: get_var_size

  integer :: i
  integer :: varid, ndims
  integer, dimension(n) :: dimids

  get_var_size(n) = -1
  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("get_var_size: can't find variable "//trim(varname))
  endif
  if (nf90_inquire_variable(ncid, varid, ndims=ndims) .ne. nf90_noerr) then
    call stop_on_err("get_var_size: can't get information for variable "//trim(varname))
  endif
  if (ndims .ne. n) then
    call stop_on_err("get_var_size: variable "//trim(varname)//" has the wrong number of dimensions")
  endif
  if (nf90_inquire_variable(ncid, varid, dimids=dimids) .ne. nf90_noerr) then
    call stop_on_err("get_var_size: can't read dimension ids for variable "//trim(varname))
  endif
  do i = 1, n
    if (nf90_inquire_dimension(ncid, dimids(i), len=get_var_size(i)) .ne. nf90_noerr) then
      call stop_on_err("get_var_size: can't get dimension lengths for variable "//trim(varname))
    endif
  enddo
end function get_var_size


!> @brief print error message and stop.
subroutine stop_on_err(msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: msg

  if (len_trim(msg) .gt. 0) then
    write(error_unit,*) trim(msg)
    stop 1
  endif
end subroutine stop_on_err


end module simple_netcdf

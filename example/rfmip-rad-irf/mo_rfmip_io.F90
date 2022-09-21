module mo_rfmip_io
use mo_rte_kind, only: wp
use mo_gas_concentrations, only: ty_gas_concs
use mo_rrtmgp_util_string, only: string_in_array
use mo_simple_netcdf, only: read_field, write_field, get_dim_size, stop_on_err
use netcdf
implicit none
private


public :: read_and_block_gases_ty
public :: read_and_block_lw_bc
public :: read_and_block_pt
public :: read_and_block_sw_bc
public :: read_size
public :: unblock_and_write


integer :: ncol_l = 0, nlay_l = 0, nexp_l = 0 ! local copies


contains


!> @brief find the size of the problem: columns, layers, perturbations (experiments)
subroutine read_size(filename, ncol, nlay, nexp)
  character(len=*), intent(in) :: filename
  integer, intent(out), optional :: ncol, nlay, nexp

  integer :: ncid

  if (nf90_open(trim(filename), nf90_nowrite, ncid) .ne. nf90_noerr) then
    call stop_on_err("read_size: can't find file "//trim(filename))
  endif

  if (present(ncol)) then
    ncol = get_dim_size(ncid, 'site')
    ncol_l = ncol
  endif
  if (present(nlay)) then
    nlay = get_dim_size(ncid, 'layer')
    nlay_l = nlay
    if (get_dim_size(ncid, 'level') .ne. nlay + 1) then
      call stop_on_err("read_size: number of levels should be nlay+1")
    endif
  endif
  if (present(nexp)) then
    nexp = get_dim_size(ncid, 'expt')
    nexp_l = nexp
  endif
  ncid = nf90_close(ncid)
end subroutine read_size


!> @brief return layer and level pressures and temperatures as arrays
!!        dimensioned (ncol, nlay/+1, nblocks)
subroutine read_and_block_pt(filename, blocksize, p_lay, p_lev, t_lay, t_lev)
  character(len=*), intent(in) :: filename
  integer, intent(in) :: blocksize
  real(kind=wp), dimension(:,:,:), allocatable, intent(out) :: p_lay, p_lev, t_lay, t_lev

  integer :: ncid
  integer :: b, nblocks
  real(kind=wp), dimension(:,:), allocatable :: temp2d
  real(kind=wp), dimension(:,:,:), allocatable :: temp3d

  if (any([ncol_l, nlay_l, nexp_l] .eq. 0)) then
    call stop_on_err("read_and_block_pt: haven't read problem size yet.")
  endif
  if (mod(ncol_l*nexp_l, blocksize) .ne. 0) then
    call stop_on_err("read_and_block_pt: number of columns doesn't fit evenly into blocks.")
  endif
  nblocks = (ncol_l*nexp_l)/blocksize
  allocate(p_lay(blocksize, nlay_l, nblocks), &
           t_lay(blocksize, nlay_l, nblocks), &
           p_lev(blocksize, nlay_l+1, nblocks), &
           t_lev(blocksize, nlay_l+1, nblocks))
  if (nf90_open(trim(filename), nf90_nowrite, ncid) .ne. nf90_noerr) then
    call stop_on_err("read_and_block_pt: can't find file "//trim(filename))
  endif

  !read p, t data; reshape to suit rrtmgp dimensions
  temp3d = reshape(spread(read_field(ncid, "pres_layer", nlay_l, ncol_l), dim=3, &
                          ncopies=nexp_l), &
                   shape=[nlay_l, blocksize, nblocks])
  do b = 1, nblocks
    p_lay(:,:,b) = transpose(temp3d(:,:,b))
  enddo
  temp3d = reshape(read_field(ncid, "temp_layer", nlay_l, ncol_l, nexp_l), &
                   shape=[nlay_l, blocksize, nblocks])
  do b = 1, nblocks
    t_lay(:,:,b) = transpose(temp3d(:,:,b))
  enddo
  deallocate(temp3d)
  temp3d = reshape(spread(read_field(ncid, "pres_level", nlay_l+1, ncol_l), dim=3, &
                          ncopies=nexp_l), &
                   shape=[nlay_l+1, blocksize, nblocks])
  do b = 1, nblocks
    p_lev(:,:,b) = transpose(temp3d(:,:,b))
  enddo
  temp3d = reshape(read_field(ncid, "temp_level", nlay_l+1, ncol_l, nexp_l), &
                   shape=[nlay_l+1, blocksize, nblocks])
  do b = 1, nblocks
    t_lev(:,:,b) = transpose(temp3d(:,:,b))
  enddo
  ncid = nf90_close(ncid)
end subroutine read_and_block_pt


!> @brief read and reshape shortwave boundary conditions.
subroutine read_and_block_sw_bc(filename, blocksize, surface_albedo, &
                                total_solar_irradiance, solar_zenith_angle)
  character(len=*), intent(in) :: filename
  integer, intent(in) :: blocksize
  real(kind=wp), dimension(:,:), allocatable, intent(out) :: surface_albedo, &
                                                             total_solar_irradiance, &
                                                             solar_zenith_angle

  integer :: ncid
  integer :: nblocks
  real(kind=wp), dimension(ncol_l, nexp_l) :: temp2d

  if (any([ncol_l, nlay_l, nexp_l] .eq. 0)) then
    call stop_on_err("read_and_block_sw_bc: haven't read problem size yet.")
  endif
  if (mod(ncol_l*nexp_l, blocksize) .ne. 0) then
    call stop_on_err("read_and_block_sw_bc: number of columns doesn't fit evenly into blocks.")
  endif
  nblocks = (ncol_l*nexp_l)/blocksize

  !check that output arrays are sized correctly : blocksize, nlay, (ncol * nexp)/blocksize
  if (nf90_open(trim(filename), nf90_nowrite, ncid) .ne. nf90_noerr) then
    call stop_on_err("read_and_block_sw_bc: can't find file "//trim(filename))
  endif
  temp2d(1:ncol_l,1:nexp_l) = spread(read_field(ncid, "surface_albedo", ncol_l), &
                                     dim=2, ncopies=nexp_l)
  surface_albedo = reshape(temp2d, shape=[blocksize, nblocks])
  temp2d(1:ncol_l,1:nexp_l) = spread(read_field(ncid, "total_solar_irradiance", ncol_l), &
                                     dim=2, ncopies=nexp_l)
  total_solar_irradiance = reshape(temp2d, shape=[blocksize, nblocks])
  temp2d(1:ncol_l,1:nexp_l) = spread(read_field(ncid, "solar_zenith_angle", ncol_l), &
                                     dim=2, ncopies=nexp_l)
  solar_zenith_angle = reshape(temp2d, shape=[blocksize, nblocks])
  ncid = nf90_close(ncid)
end subroutine read_and_block_sw_bc


!> @brief  read and reshape longwave boundary conditions.
subroutine read_and_block_lw_bc(filename, blocksize, surface_emissivity, surface_temperature)
  character(len=*), intent(in) :: filename
  integer, intent(in) :: blocksize
  real(kind=wp), dimension(:,:), allocatable, intent(out) :: surface_emissivity
  real(kind=wp), dimension(:,:), allocatable, intent(out) :: surface_temperature

  integer :: ncid
  integer :: nblocks
  real(kind=wp), dimension(ncol_l, nexp_l) :: temp2d

  if (any([ncol_l, nlay_l, nexp_l] .eq. 0)) then
    call stop_on_err("read_and_block_lw_bc: haven't read problem size yet.")
  endif
  if (mod(ncol_l*nexp_l, blocksize) .ne. 0) then
    call stop_on_err("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
  endif
  nblocks = (ncol_l*nexp_l)/blocksize

  if (nf90_open(trim(filename), nf90_nowrite, ncid) .ne. nf90_noerr) then
    call stop_on_err("read_and_block_lw_bc: can't find file "//trim(filename))
  endif

  !allocate on assigment
  temp2d(1:ncol_l,1:nexp_l) = spread(read_field(ncid, "surface_emissivity", ncol_l), &
                                     dim=2, ncopies=nexp_l)
  surface_emissivity = reshape(temp2d, shape=[blocksize, nblocks])
  temp2d(1:ncol_l,1:nexp_l) = read_field(ncid, "surface_temperature", ncol_l, nexp_l)
  surface_temperature = reshape(temp2d, shape=[blocksize, nblocks])
  ncid = nf90_close(ncid)
end subroutine read_and_block_lw_bc


!> @brief read and reshape gas concentrations.
subroutine read_and_block_gases_ty(filename, blocksize, gas_names, names_in_file, &
                                   gas_conc_array)
  character(len=*), intent(in) :: filename
  integer, intent(in) :: blocksize
  character(len=*), dimension(:), intent(in) :: gas_names
  character(len=*), dimension(:), intent(in) :: names_in_file
  type(ty_gas_concs), dimension(:), allocatable, intent(out) :: gas_conc_array

  integer :: ncid
  integer :: nblocks
  integer :: b, g
  integer, dimension(:,:), allocatable :: exp_num
  real(kind=wp), dimension(:), allocatable :: gas_conc_temp_1d
  real(kind=wp), dimension(:,:,:), allocatable :: gas_conc_temp_3d

  if (any([ncol_l, nlay_l, nexp_l] .eq. 0)) then
    call stop_on_err("read_and_block_lw_bc: haven't read problem size yet.")
  endif
  if (mod(ncol_l*nexp_l, blocksize) .ne. 0) then
    call stop_on_err("read_and_block_lw_bc: number of columns doesn't fit evenly into blocks.")
  endif
  nblocks = (ncol_l*nexp_l)/blocksize
  allocate(gas_conc_array(nblocks))
  do b = 1, nblocks
    if (any([(string_in_array(gas_names(g), ['h2o', 'o3 ', 'no2']), g=1, size(gas_names))])) then
      call stop_on_err(gas_conc_array(b)%init(gas_names))
    else
      call stop_on_err(gas_conc_array(b)%init([gas_names, 'h2o    ', 'o3     ', 'no2    ']))
    endif
  enddo

  !experiment index for each colum
  exp_num = reshape(spread([(b, b = 1, nexp_l)], 1, ncopies=ncol_l), &
                    shape=[blocksize, nblocks], order=[1,2])

  if (nf90_open(trim(filename), nf90_nowrite, ncid) .ne. nf90_noerr) then
    call stop_on_err("read_and_block_gases_ty: can't find file "//trim(filename))
  endif

  !water vapor and ozone depend on col, lay, exp: look just like other fields
  gas_conc_temp_3d = reshape(read_field(ncid, "water_vapor", nlay_l, ncol_l, nexp_l), &
                             shape=[nlay_l, blocksize, nblocks])* &
                     read_scaling(ncid, "water_vapor")
  do b = 1, nblocks
    call stop_on_err(gas_conc_array(b)%set_vmr('h2o', transpose(gas_conc_temp_3d(:,:,b))))
  enddo

  gas_conc_temp_3d = reshape(read_field(ncid, "ozone", nlay_l, ncol_l, nexp_l), &
                             shape=[nlay_l, blocksize, nblocks])* &
                     read_scaling(ncid, "ozone")
  do b = 1, nblocks
    call stop_on_err(gas_conc_array(b)%set_vmr('o3', transpose(gas_conc_temp_3d(:,:,b))))
  enddo

  !all other gases are a function of experiment only
  do g = 1, size(gas_names)
    !skip 3d fields above, also no2 since rfmip doesn't have this
    if (string_in_array(gas_names(g), ['h2o', 'o3 ', 'no2'])) then
      cycle
    endif

    !read the values as a function of experiment
    gas_conc_temp_1d = read_field(ncid, trim(names_in_file(g))//"_GM", nexp_l)* &
                       read_scaling(ncid, trim(names_in_file(g))//"_GM")

    do b = 1, nblocks
      !does every value in this block belong to the same experiment?
      if (all(exp_num(1,b) .eq. exp_num(2:,b))) then
        !provide a scalar value
        call stop_on_err(gas_conc_array(b)%set_vmr(gas_names(g), &
                                                   gas_conc_temp_1d(exp_num(1,b))))
      else
        !create 2d field, blocksize x nlay, with scalar values from each experiment
        call stop_on_err(gas_conc_array(b)%set_vmr(gas_names(g), &
                                                   spread(gas_conc_temp_1d(exp_num(:,b)), 2, &
                                                          ncopies=nlay_l)))
      endif
    enddo

    !no2 is the one gas known to the k-distribution that isn't provided by rfmip
    !it would be better to remove it from
    do b = 1, nblocks
      call stop_on_err(gas_conc_array(b)%set_vmr('no2', 0._wp))
    enddo
  enddo
  ncid = nf90_close(ncid)
end subroutine read_and_block_gases_ty


function read_scaling(ncid, varname)
  integer, intent(in) :: ncid
  character(len=*), intent(in) :: varname
  real(kind=wp) :: read_scaling

  integer :: varid
  character(len=16) :: charunits

  if (nf90_inq_varid(ncid, trim(varname), varid) .ne. nf90_noerr) then
    call stop_on_err("read_scaling: can't find variable "//trim(varname))
  endif
  if (nf90_get_att(ncid, varid, "units", charunits) .ne. nf90_noerr) then
    call stop_on_err("read_scaling: can't read attribute 'units' from variable "//trim(varname))
  endif
  read(charunits, *) read_scaling
  return
end function read_scaling


!> @brief reshape and reorder values (nominally fluxes) from rte order (ncol, nlev, nblocks)
!!        to rfmip order (nlev, ncol, nexp), then write them to a user-specified variable
!!        in a netcdf file.
subroutine unblock_and_write(filename, varname, values)
  character(len=*), intent(in) :: filename, varname
  real(kind=wp), dimension(:,:,:), intent(in) :: values

  integer :: ncid
  integer :: b, blocksize, nlev, nblocks
  real(kind=wp), dimension(:,:), allocatable :: temp2d

  if (any([ncol_l, nlay_l, nexp_l] .eq. 0)) then
    call stop_on_err("unblock_and_write: haven't read problem size yet.")
  endif
  blocksize = size(values, 1)
  nlev = size(values, 2)
  nblocks = size(values, 3)
  if (nlev .ne. nlay_l + 1) then
    call stop_on_err('unblock_and_write: array values has the wrong number of levels')
  endif
  if (blocksize*nblocks .ne. ncol_l*nexp_l) then
    call stop_on_err('unblock_and_write: array values has the wrong number of blocks/size')
  endif
  allocate(temp2d(nlev, ncol_l*nexp_l))
  do b = 1, nblocks
    temp2d(1:nlev, ((b-1)*blocksize+1):(b*blocksize)) = transpose(values(1:blocksize,1:nlev,b))
  enddo
  if (nf90_open(trim(filename), nf90_write, ncid) .ne. nf90_noerr) then
    call stop_on_err("unblock_and_write: can't find file "//trim(filename))
  endif
  call stop_on_err(write_field(ncid, varname, reshape(temp2d, shape=[nlev, ncol_l, nexp_l])))
  ncid = nf90_close(ncid)
end subroutine unblock_and_write


end module mo_rfmip_io

!--------------------------------------------------------------------------------------------------------------------
module mo_testing_io
  !
  ! RTE+RRTMGP modules
  !
  use mo_rte_kind,           only: wp
  !
  ! NetCDF I/O routines, shared with other RTE+RRTMGP examples
  !
  use mo_simple_netcdf,      only: write_field, create_dim, create_var
  use netcdf
  implicit none
  private
  public :: write_broadband_field
contains
  !--------------------------------------------------------------------------------------------------------------------
  subroutine write_broadband_field(fileName, field, field_name, field_description, col_dim_name, vert_dim_name)
    !
    ! Write a field defined by column and some vertical dimension (lev or lay))
    !
    character(len=*),         intent(in) :: fileName
    real(wp), dimension(:,:), intent(in) :: field
    character(len=*),         intent(in) :: field_name, field_description
    character(len=*), optional, &
                              intent(in) ::col_dim_name, vert_dim_name
    ! -------------------
    integer           :: ncid, varid, ncol, nlev
    !
    ! Names of column (first) and vertical (second) dimension.
    !   Because they are used in an array constuctor the need to have the same number of characters
    !
    character(len=32) :: cdim, vdim
    ! -------------------
    cdim = "site "
    vdim = "level"
    if(present(col_dim_name))  cdim = col_dim_name
    if(present(vert_dim_name)) vdim = vert_dim_name

    if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

    ncol  = size(field, dim=1)
    nlev  = size(field, dim=2)

    call create_dim(ncid, trim(cdim), ncol)
    call create_var(ncid, trim(field_name),  [cdim, vdim], [ncol, nlev])
    call stop_on_err(write_field(ncid, trim(field_name),  field))
    !
    ! Adding descriptive text as an attribute means knowing the varid
    !
    if(nf90_inq_varid(ncid, trim(field_name), varid) /= NF90_NOERR) &
      call stop_on_err("Can't find variable " // trim(field_name))
    if(nf90_put_att(ncid, varid, "description", trim(field_description)) /= NF90_NOERR) &
      call stop_on_err("Can't write 'description' attribute to variable " // trim(field_name))

    ncid = nf90_close(ncid)

  end subroutine write_broadband_field
  !--------------------------------------------------------------------------------------------------------------------
end module mo_testing_io

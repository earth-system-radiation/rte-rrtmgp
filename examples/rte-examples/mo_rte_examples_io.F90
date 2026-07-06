module mo_rte_examples_io
  use mo_rte_kind,           only: wp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_gas_optics_util_string, &
                            only: lower_case, string_in_array
  use mo_testing_utils,     only: stop_on_err
  use mo_simple_netcdf,     only: get_dim_size, read_field
  use netcdf
  implicit none

  private
  public :: inquire_rte_example, read_rte_example
  ! h2o, co2, ch4, n2o, co, n2, o2, cfc11, cfc12, ! DDQ
  ! h2o co2 o3 n2o co ch4 o2 n2 ccl4 cfc11 cfc12 cfc22 hfc143a hfc125 hfc23 hfc32 hfc134a cf4 no2 ! all gases in RRTMGP
  character(len=24), parameter, &
    dimension(9) :: state_variables = [ &
    "pres_layer            ", &
    "pres_level            ", &
    "temp_layer            ", &
    "temp_level            ", &
    "surface_albedo        ", &
    "solar_zenith_angle    ", &
    "total_solar_irradiance", &
    "surface_emissivity    ", &
    "surface_temperature   "  ]

  ! Module variable - set when data is read, used when writing out resuls
  integer :: nvar = 0
contains
  ! --------------------------------------------------
  !
  ! Get problem sizes and the list of available gas concentrations (which RRTMGP needs)
  ! Anything that isn't a state variable is assumed to be a gas concentration
  ! State variables are
  !    pres_layer, pres_leve, temp_layer, temp_level, ps?
  !    surface_albedo, solar_zenith_angle, total_solar_irradiance
  !    surface_emissivity, surface_temperature
  !
  subroutine inquire_rte_example(filename, &
  	               ncol, nlay, available_gases)
    character(len=*),   intent(in)  :: filename
    integer,            intent(out) :: ncol, nlay
    type(ty_gas_concs), intent(out) :: available_gases

    integer :: nvariant = 1
    integer :: ncid, ivar, varid
    character(len=32) :: varname
    ! -----
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("inquire_rte_example: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'layer')
    nvar = get_dim_size(ncid, 'variant')

    call available_gases%reset()
    ivar = 1
    do
      if (nf90_inquire_variable(ncid, ivar, name = varname) /= NF90_NOERR) then
        exit
      else
        if (.not. string_in_array(trim(varname), state_variables)) &
          call stop_on_err(available_gases%set_vmr(trim(varname), 0._wp))
      end if
      ivar = ivar + 1
    end do

    ncid = nf90_close(ncid)

  end subroutine inquire_rte_example
  ! --------------------------------------------------
  subroutine read_rte_example(filename, is_lw, &
  	               pres_layer, pres_level, temp_layer, temp_level, gas_concs, &
  	               surface_albedo, solar_zenith_angle, total_solar_irradiance, &
  	               surface_emissivity, surface_temperature)
    character(len=*), intent(in) :: filename
    logical,          intent(in) :: is_lw
    real(wp), intent(out), dimension(:,:), allocatable &
       :: pres_layer, pres_level, temp_layer, temp_level
    type(ty_gas_concs), intent(out) :: gas_concs
    real(wp), intent(out), dimension(:), allocatable , &
      optional :: surface_albedo, solar_zenith_angle, total_solar_irradiance, &
                   surface_emissivity, surface_temperature

    integer :: ncid
    integer :: ncol, nlay, nvar, ngas, igas
    character(len=32), dimension(:), allocatable :: gas_names
    ! -----
    if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_rte_example: can't find file " // trim(fileName))

    ncol = get_dim_size(ncid, 'col')
    nlay = get_dim_size(ncid, 'layer')
    nvar = max(1, get_dim_size(ncid, 'variant'))
    !
    ! State variables
    !
    pres_layer = reorder(read_field(ncid, "pres_layer", &
    	                 ncol, nlay,   nvar))
    pres_level = reorder(read_field(ncid, "pres_level", &
    	                 ncol, nlay+1, nvar))
    temp_layer = reorder(read_field(ncid, "temp_layer", &
    	                 ncol, nlay,   nvar))
    temp_level = reorder(read_field(ncid, "temp_level", &
    	                 ncol, nlay+1, nvar))
    if(present(surface_emissivity)) &
	    surface_emissivity = reshape(read_field(ncid, "surface_emissivity", &
	    	                 ncol, nvar), &
	    	                 shape = [ncol * nvar])
    if(present(surface_temperature)) &
	    surface_temperature = reshape(read_field(ncid, "surface_temperature", &
	    	                 ncol, nvar), &
	    	                 shape = [ncol * nvar])
    if(present(surface_albedo)) &
	    surface_albedo = reshape(read_field(ncid, "surface_albedo", &
	    	                 ncol, nvar), &
	    	                 shape = [ncol * nvar])
    if(present(solar_zenith_angle)) &
	    solar_zenith_angle = reshape(read_field(ncid, "solar_zenith_angle", &
	    	                 ncol, nvar), &
	    	                 shape = [ncol * nvar])
    if(present(total_solar_irradiance)) &
	    total_solar_irradiance = reshape(read_field(ncid, "total_solar_irradiance", &
	    	                 ncol, nvar), &
	    	                 shape = [ncol * nvar])

    gas_names = gas_concs%get_gas_names()
    do igas = 1, gas_concs%get_num_gases()
          call stop_on_err(                    &
          	gas_concs%set_vmr(gas_names(igas), &
          	                  reorder(read_field(ncid, gas_names(igas), &
    	                              ncol, nlay, nvar)) )              &
          )
    end do

    ncid = nf90_close(ncid)

  end subroutine read_rte_example
  ! --------------------------------------------------
  !
  ! Convert from ncol, nlay, nvariant to ncol*nvariant, nlay
  !
  pure function reorder(array_in) result(array_out)
    real(wp), dimension(:,:,:), intent(in) :: array_in
    real(wp), dimension(size(array_in, 1) * size(array_in, 3), &
    	                size(array_in, 2)) :: array_out

    integer :: nx, ny, nv, v

    nx = size(array_in, 1)
    ny = size(array_in, 2)
    nv = size(array_in, 3)

    do v = 1, nv
      array_out(nx*(v-1)+1:nx*v, :) = array_in(:, :, v)
    end do
  end function reorder
  ! --------------------------------------------------
end module mo_rte_examples_io

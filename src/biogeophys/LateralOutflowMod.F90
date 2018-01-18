module LateralOutflowMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Type and associated routines for calculating lateral outflow
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use abortutils   , only : endrun
  use clm_varctl   , only : iulog
  use clm_varcon   , only : e_ice, rpi
  use clm_varpar   , only : nlevsoi
  use ColumnType   , only : column_type
  use GridcellType , only : gridcell_type
  use SoilHydrologyType, only : soilhydrology_type

  implicit none
  save
  private

  ! !PRIVATE DATA MEMBERS:

  ! Indices into baseflow_methods array
  integer, parameter :: METHOD_INDEX_SOIL  = 1
  integer, parameter :: METHOD_INDEX_CROP  = 2
  integer, parameter :: METHOD_INDEX_OTHER = 3
  integer, parameter :: NUM_METHODS = 3

  ! !PUBLIC TYPES:

  type, public :: lateral_outflow_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: qflx_latflow_out_col(:)     ! lateral flow output (mm/s)
     real(r8), pointer, public :: qflx_latflow_out_vol_col(:) ! lateral flow output volume (m^3/s)

     ! Private data members
     integer :: baseflow_methods(NUM_METHODS) ! baseflow method to use over various landunits

     integer :: transmissivity_method  ! Only used for BASEFLOW_METHOD_KINEMATIC and BASEFLOW_METHOD_DARCY

     real(r8) :: baseflow_scalar  ! Only used for BASEFLOW_METHOD_POWER_LAW
   contains
     ! Public routines
     procedure, public :: Init
     procedure, public :: LateralOutflow ! Calculate lateral outflow

     ! Private routines
     procedure, private :: ReadNamelist
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

     procedure, private :: ComputeLateralOutflowPowerLaw
     procedure, private, nopass :: ConvertLatflowOutToVolume      ! simple method to convert a latflow output flux to a volume
  end type lateral_outflow_type

  ! !PRIVATE DATA MEMBERS:

  integer, parameter :: BASEFLOW_METHOD_POWER_LAW = 1

  integer, parameter :: TRANSMISSIVITY_METHOD_LAYERSUM = 1
  integer, parameter :: TRANSMISSIVITY_METHOD_CONSTANT = 2
  integer, parameter :: TRANSMISSIVITY_METHOD_POWER    = 3

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  ! ========================================================================
  ! Infrastructure routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds, NLFilename)
    !
    ! !DESCRIPTION:
    ! Initialize this lateral_outflow_type object
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    character(len=*), intent(in) :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------

    call this%ReadNamelist(NLFilename)
    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine ReadNamelist(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the lateraloutflow namelist
    !
    ! !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(inout) :: this
    character(len=*), intent(in) :: NLFilename ! namelist filename
    !
    ! !LOCAL VARIABLES:

    ! temporary variables corresponding to the variables read from namelist
    real(r8) :: baseflow_scalar

    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=*), parameter :: nmlname = 'lateraloutflow_inparm'

    character(len=*), parameter :: subname = 'ReadNamelist'
    !-----------------------------------------------------------------------

    namelist /lateraloutflow_inparm/ baseflow_scalar

    ! Initialize parameters to garbage defaults, forcing all to be specified explicitly in
    ! order to get reasonable results
    baseflow_scalar = nan

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=lateraloutflow_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (baseflow_scalar, mpicom)

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=lateraloutflow_inparm)
       write(iulog,*) ' '
    end if

    this%baseflow_scalar = baseflow_scalar

    ! TODO(wjs, 2018-01-05) We'll read these from namelist. (Probably read three separate
    ! namelist variables: baseflow_method_soil, baseflow_method_crop,
    ! baseflow_method_other.)
    this%baseflow_methods(METHOD_INDEX_SOIL)  = BASEFLOW_METHOD_POWER_LAW
    this%baseflow_methods(METHOD_INDEX_CROP)  = BASEFLOW_METHOD_POWER_LAW
    this%baseflow_methods(METHOD_INDEX_OTHER) = BASEFLOW_METHOD_POWER_LAW

    this%transmissivity_method = TRANSMISSIVITY_METHOD_LAYERSUM

  end subroutine ReadNamelist


  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Allocate memory for this lateral_outflow_type object
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate(this%qflx_latflow_out_col(begc:endc))     ; this%qflx_latflow_out_col(:)     = nan
    allocate(this%qflx_latflow_out_vol_col(begc:endc)) ; this%qflx_latflow_out_vol_col(:) = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize lateral_outflow_type history variables
    !
    ! !USES:
    use histFileMod , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    ! Nothing to do for now

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Perform cold-start initialization for lateral_outflow_type
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    ! Nothing to do

  end subroutine InitCold

  ! ========================================================================
  ! Science routines
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine LateralOutflow(this, bounds, num_hydrologyc, filter_hydrologyc, &
       col, grc, soilhydrology_inst, jwt, dzmm)
    !
    ! !DESCRIPTION:
    ! Compute lateral outflow
    !
    ! This has been described elsewhere as 'topographic runoff'
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(inout) :: this
    type(bounds_type)        , intent(in) :: bounds               
    integer                  , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer                  , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(column_type)        , intent(in) :: col
    type(gridcell_type)      , intent(in) :: grc
    type(soilhydrology_type) , intent(in) :: soilhydrology_inst
    integer                  , intent(in) :: jwt( bounds%begc: )  ! index of the soil layer right above the water table (-)
    real(r8)                 , intent(in) :: dzmm( bounds%begc: , 1: ) ! layer thickness (mm)
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'LateralOutflow'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(jwt) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dzmm) == (/bounds%endc, nlevsoi/)), errMsg(sourcefile, __LINE__))

    associate( &
         qflx_latflow_out     => this%qflx_latflow_out_col     , & ! Output: [real(r8) (:) ] lateral flow output (mm/s)
         qflx_latflow_out_vol => this%qflx_latflow_out_vol_col   & ! Output: [real(r8) (:) ] lateral flow output volume (m^3/s)
         )

    ! FIXME(wjs, 2018-01-03) divide the filter into two sub-filters based on whether
    ! istsoil is true
    !
    ! To do this, introduce a new routine into filterColMod: col_filter_divide_using_ltype.
    ! This will be a cross between col_filter_from_ltypes and
    ! col_filter_from_filter_and_logical_array; it will accept an existing filter and an
    ! integer giving ltype; it will return TWO filters whose union is the original filter.
    !
    ! Update (based on input from the ctsm meeting 2018-01-04): Actually, do three: one
    ! for istsoil, one for istcrop, and one for everything else. To do this, I think I
    ! should have a routine that divides an existing filter by landunit types: its inputs
    ! are an existing filter and an array of integers giving ltypes; its output is an
    ! array of filters (n+1 filters, where n is the number of ltypes): filter 1 applies
    ! over ltype 1, filter 2 applies over ltype 2... filter n+1 applies over
    ! non-specified ltypes.
    !
    !   Then in this routine, I'll have an array of filters. I'll have parameters saying
    !   which index of this array of filters is what. e.g., filter_istsoil = 1,
    !   filter_istcrop = 2, filter_other = 3.
    !
    !   The implementation of the filter-creator can have, at the beginning, a creation
    !   of an array mapping ltypes to the appropriate index in the output array. e.g.,
    !
    !     do i = 1, max_lunit
    !       filter_num(i) = [whatever code is needed to determine the right filter number
    !       for this landunit]
    !     end do
    !
    !     Then loop over the original filter, assigning points to the correct output
    !     filter.

    ! FIXME(wjs, 2018-01-04) Then have a select case for each sub-filter, with the
    ! appropriate namelist item.

    call this%ComputeLateralOutflowPowerLaw(bounds, num_hydrologyc, filter_hydrologyc, &
         col, grc, soilhydrology_inst, &
         jwt = jwt(bounds%begc:bounds%endc), &
         dzmm = dzmm(bounds%begc:bounds%endc,:), &
         baseflow_scalar = this%baseflow_scalar, &
         qflx_latflow_out = qflx_latflow_out(bounds%begc:bounds%endc), &
         qflx_latflow_out_vol = qflx_latflow_out_vol(bounds%begc:bounds%endc))

    end associate

  end subroutine LateralOutflow

  !-----------------------------------------------------------------------
  subroutine ComputeLateralOutflowPowerLaw(this, bounds, num_c, filter_c, &
       col, grc, soilhydrology_inst, jwt, dzmm, baseflow_scalar, &
       qflx_latflow_out, qflx_latflow_out_vol)
    !
    ! !DESCRIPTION:
    ! Compute lateral outflow using a power law method
    !
    ! !ARGUMENTS:
    class(lateral_outflow_type), intent(in) :: this
    type(bounds_type)        , intent(in) :: bounds               
    integer                  , intent(in) :: num_c       ! number of column points in column filter
    integer                  , intent(in) :: filter_c(:) ! column filter
    type(column_type)        , intent(in) :: col
    type(gridcell_type)      , intent(in) :: grc
    type(soilhydrology_type) , intent(in) :: soilhydrology_inst
    integer                  , intent(in) :: jwt( bounds%begc: )  ! index of the soil layer right above the water table (-)
    real(r8)                 , intent(in) :: dzmm( bounds%begc: , 1: ) ! layer thickness (mm)
    real(r8)                 , intent(in) :: baseflow_scalar

    ! The following are set over the given filter, and are left unchanged elsewhere
    real(r8) , intent(inout) :: qflx_latflow_out( bounds%begc: )     ! lateral flow output (mm/s)
    real(r8) , intent(inout) :: qflx_latflow_out_vol( bounds%begc: ) ! lateral flow output volume (m^3/s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc ! column filter index
    integer :: c  ! column index
    integer :: j  ! level index

    real(r8) :: dzsum           ! summation of dzmm of layers below water table (mm)
    real(r8) :: icefracsum      ! summation of icefrac*dzmm of layers below water table (-)
    real(r8) :: imped

    real(r8), parameter :: n_baseflow = 1 !drainage power law exponent

    character(len=*), parameter :: subname = 'ComputeLateralOutflowPowerLaw'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(jwt) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(dzmm) == (/bounds%endc, nlevsoi/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_latflow_out) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_latflow_out_vol) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    associate( &
         zi         =>    col%zi                         , & ! Input: [real(r8) (:,:) ] interface level below a "z" level (m)           
         nbedrock   =>    col%nbedrock                   , & ! Input: [real(r8) (:,:) ] depth to bedrock (m)           
         topo_slope =>    col%topo_slope                 , & ! Input: [real(r8) (:)   ] topographic slope
         icefrac    =>    soilhydrology_inst%icefrac_col , & ! Input: [real(r8) (:,:) ] fraction of ice in layer
         zwt        =>    soilhydrology_inst%zwt_col       & ! Input: [real(r8) (:)   ] water table depth (m)                             
         )

    do fc = 1, num_c
       c = filter_c(fc)

       dzsum = 0._r8
       icefracsum = 0._r8
       do j = max(jwt(c),1), nlevsoi
          dzsum  = dzsum + dzmm(c,j)
          icefracsum = icefracsum + icefrac(c,j) * dzmm(c,j)
       end do
       imped=10._r8**(-e_ice*(icefracsum/dzsum))
       !@@
       ! baseflow is power law expression relative to bedrock layer
       if(zwt(c) <= zi(c,nbedrock(c))) then 
          qflx_latflow_out(c) = imped * baseflow_scalar * tan(rpi/180._r8*topo_slope(c))* &
               (zi(c,nbedrock(c)) - zwt(c))**(n_baseflow)
          if (qflx_latflow_out(c) < 0._r8) then
             write(iulog,*) subname//' ERROR: Unexpected negative qflx_latflow_out'
             write(iulog,*) 'c, qflx_latflow_out = ', c, qflx_latflow_out(c)
             call endrun(msg=subname//' ERROR: Unexpected negative qflx_latflow_out ' // &
                  errMsg(sourcefile, __LINE__))
          end if
       else
          qflx_latflow_out(c) = 0._r8
       endif

    end do

    call this%ConvertLatflowOutToVolume(bounds, num_c, filter_c, &
         col, grc, &
         qflx_latflow_out = qflx_latflow_out(bounds%begc:bounds%endc), &
         qflx_latflow_out_vol = qflx_latflow_out_vol(bounds%begc:bounds%endc))

    end associate

  end subroutine ComputeLateralOutflowPowerLaw

  !-----------------------------------------------------------------------
  subroutine ConvertLatflowOutToVolume(bounds, num_c, filter_c, &
       col, grc, qflx_latflow_out, qflx_latflow_out_vol)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in) :: bounds               
    integer                  , intent(in) :: num_c       ! number of column points in column filter
    integer                  , intent(in) :: filter_c(:) ! column filter
    type(column_type)        , intent(in) :: col
    type(gridcell_type)      , intent(in) :: grc

    real(r8) , intent(in)    :: qflx_latflow_out( bounds%begc: )     ! lateral flow output (mm/s)
    ! The following is set over the given filter, and is left unchanged elsewhere
    real(r8) , intent(inout) :: qflx_latflow_out_vol( bounds%begc: ) ! lateral flow output volume (m^3/s)
    !
    ! !LOCAL VARIABLES:
    integer :: fc ! column filter index
    integer :: c  ! column index
    integer :: g  ! gridcell index

    real(r8), parameter :: mm_to_m = 1.e-3
    real(r8), parameter :: km2_to_m2 = 1.e6

    character(len=*), parameter :: subname = 'ConvertLatflowOutToVolume'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(qflx_latflow_out) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    SHR_ASSERT_ALL((ubound(qflx_latflow_out_vol) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))

    do fc = 1, num_c
       c = filter_c(fc)
       g = col%gridcell(c)
       qflx_latflow_out_vol(c) = qflx_latflow_out(c) * mm_to_m * &
            (col%wtgcell(c) * grc%area(g) * km2_to_m2)
    end do


  end subroutine ConvertLatflowOutToVolume


end module LateralOutflowMod

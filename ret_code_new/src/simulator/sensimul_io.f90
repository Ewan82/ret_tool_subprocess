!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file sensimul_io.f90
!> \brief provides I/O handling for combined S1/S2 retrieval system
!> \authors The Inversion Lab (Michael Vossbeck)
!> \date  April, October 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!***********************************************************
!     ncrd_retrieval_config
!
!> @brief Prepares required basic setup for running the Sentinel Simulator observational
!         operator. These informations must be provided by a NetCDF file.
!
!> @details TBD
!
!> @param[in]  fname        file name of NetCDF files providing setup information
!> @param[out] succeed      flag indicating whether setup was successfully read
!
subroutine ncrd_retrieval_config(fname, succeed)
  use netcdf
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname
  logical, intent(out) :: succeed
  ! local decls
  integer :: fid,rc

  !-- initialise output flag
  succeed = .false.

  !-- open file handle
  rc = nf90_open(trim(fname), nf90_nowrite, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config:opening file failed',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- get 'npoints' handle
  if( .not. set_npoints() ) then
     return
  endif

  !-- simulation type (S1,S2)
  if( .not. load_timept_simtyp('sim_typ') ) then
     return
  endif

  !-- determein S1,S2 indices

  !-- illumination-view geometry (for S2 only)
  if( .not. load_ivgeom('ivgeom') ) then
     return
  endif

  if( .not. load_timevalues() ) then
     return
  endif

  !-- eventually successfull load of configuration
  succeed = .true.

  !-- close handle
  rc = nf90_close(fid)

  contains
    logical function set_npoints()
      implicit none
      ! local decls
      integer :: dimid,dimval,rcc

      set_npoints = .false.

      !-- get 'npoints' handle
      rcc = nf90_inq_dimid(fid, 'npoints', dimid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config:'//&
              'npoints dimension not accessible',&
              trim(nf90_strerror(rcc))
         return
      else
         rcc = nf90_inquire_dimension(fid, dimid, len=dimval)
         if( rcc.ne.nf90_noerr ) then
            write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config:'//&
                 'npoints dimension extent',&
                 trim(nf90_strerror(rcc))
            return
         else
            npts = dimval
         endif
      endif

      set_npoints = .true.
    end function set_npoints

    logical function load_timept_simtyp(dsname)
      implicit none
      ! arguments
      character(len=*), intent(in) :: dsname
      ! local decls
      integer :: vid, rcc
      integer :: i,i_s1,i_s2

      !-- initialise return value
      load_timept_simtyp = .false.

      !-- allocate global structure
      allocate( timept_simtyp(npts) )
      rcc = nf90_inq_varid(fid, dsname, vid)
      if( rcc.ne. nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::',trim(nf90_strerror(rcc))
         return
      end if

      !-- read time_typ data
      rcc = nf90_get_var(fid, vid, timept_simtyp)
      if( rcc.ne. nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::',trim(nf90_strerror(rcc))
         return
      end if

      !-- determine mask for S1 simulations
      npts_s1 = count( iand(timept_simtyp,1).eq.1 )

      !-- determine mask for S2 simulations
      npts_s2 = count( iand(timept_simtyp,2).eq.2 )

      !-- setup index arrays for S1,S2
      if( allocated(timept_idxs_s1) ) then
         write(*, '(a)') ' FATAL::load_timept_simtyp:'//&
              'array timept_idxs_s1 already allocated when entering initial setup!'
         return         
      else
         allocate(timept_idxs_s1(npts_s1))
      endif
      if( allocated(timept_idxs_s2) ) then
         write(*, '(a)') ' FATAL::load_timept_simtyp:'//&
              'array timept_idxs_s2 already allocated when entering initial setup!'
         return         
      else
         allocate(timept_idxs_s2(npts_s2))
      endif
      i_s1 = 1
      i_s2 = 1
      do i=1,npts
         if( idx_is_s1(i) ) then
            timept_idxs_s1(i_s1) = i
            i_s1 = i_s1+1
         endif
         if( idx_is_s2(i) ) then
            timept_idxs_s2(i_s2) = i
            i_s2 = i_s2+1
         endif
      enddo

      !-- eventually successful return
      load_timept_simtyp = .true.

    end function load_timept_simtyp

    logical function load_ivgeom(dsname)
      implicit none
      ! arguments
      character(len=*), intent(in) :: dsname
      ! local decls
      integer :: vid, rcc
      integer :: ndims, dimids(nf90_max_var_dims)
      integer :: dimsiz
      logical :: ldebug = .true.

      !-- initialise return value
      load_ivgeom = .false.

      rcc = nf90_inq_varid(fid, dsname, vid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::load_ivgeom:',trim(nf90_strerror(rcc))
         return
      end if
      rcc = nf90_inquire_variable(fid, vid, ndims=ndims, dimids=dimids)
      if( rcc.ne.nf90_noerr ) then
         write(*,'(a)') ' ###NETCDF-ERROR::load_ivgeom:'//&
              'variable dimensions could not be accessed'
         return
      else if( ndims.ne.2 ) then
         write(*, '(a,i2)') ' FATAL::load_ivgeom:'//&
              'expected 2-dimensional array, but got ndims=',ndims
         return
      else
         if(ldebug) then
            write(*,'(a,i2)') ' DEBUG::load_ivgeom:ndims=',ndims
         endif
      endif
      rcc = nf90_inquire_dimension(fid, dimids(1), len=dimsiz)
      if( dimsiz.ne.4 ) then
         write(*, '(a,i3)') ' FATAL::load_ivgeom:'//&
              'first dimension size expected=4, got=',dimsiz
         return
      endif
      rcc = nf90_inquire_dimension(fid, dimids(2), len=dimsiz)
      if( dimsiz.ne.npts ) then
         write(*, '(a,2(a,i3,1x))') ' FATAL::load_ivgeom:'//&
              'second dimension size ', 'expected=',npts,'got=',dimsiz
      endif

      if( allocated(iv_geom) ) then
         write(*, '(a)') ' FATAL::load_ivgeom:'//&
              'array iv_geom already allocated when entering initial setup!'
         return
      else
         allocate( iv_geom(4,npts) )
      endif
      !-- read illumination-view data
      rcc = nf90_get_var(fid, vid, iv_geom)
      if( rcc.ne.nf90_noerr ) then
         write(*,'(a)') ' NETCDF-ERROR::load_ivgeom:reading iv_geom failed. '//&
              ' ('//trim(nf90_strerror(rcc))//')'
         return
      endif

      !-- successfully loaded
      load_ivgeom = .true.

    end function load_ivgeom


    logical function load_timevalues()
      implicit none
      ! local decls
      character(len=*), parameter :: dsname = 'time'
      integer :: vid, rcc
      integer :: ndims, dimids(nf90_max_var_dims)
      integer :: dimsiz, attr_len
      logical :: ldebug = .false.


      load_timevalues = .false.

      rcc = nf90_inq_varid(fid, dsname, vid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::load_timevalues:',trim(nf90_strerror(rcc))
         return
      end if
      rcc = nf90_inquire_variable(fid, vid, ndims=ndims, dimids=dimids)
      if( rcc.ne.nf90_noerr ) then
         write(*,'(a)') ' ###NETCDF-ERROR::load_timevalues:'//&
              'variable dimensions could not be accessed'
         return
      else if( ndims.ne.1 ) then
         write(*, '(a,i2)') ' FATAL::load_timevalues:'//&
              'expected 1-dimensional array, but got ndims=',ndims
         return
      else
         if(ldebug) then
            write(*,'(a,i2)') ' DEBUG::load_timevalues:ndims=',ndims
         endif
      endif
      rcc = nf90_inquire_dimension(fid, dimids(1), len=dimsiz)
      if( dimsiz.ne.npts ) then
         write(*, '(a,2(a,i4,1x))') ' FATAL::load_timevalues:'//&
              'unexpected dimension size:', 'expected=',npts,'got=',dimsiz
         return
      endif

      !-- allocate global time variable
      if( allocated(timepts) ) then
         write(*, '(a)') ' FATAL::load_timevalues:'//&
              'array ---timepts--- already allocated when entering initial I/O'
         return
      else
         allocate( timepts(npts) )
      endif

      !-- read time values
      rcc = nf90_get_var(fid, vid, timepts)
      if( rcc.ne.nf90_noerr ) then
         write(*, '(a)') ' ###NETCDF-ERROR::load_timevalues:'//&
              'error while reading dataset.'
         return
      endif

      !-- determine the unit
      rcc = nf90_inquire_attribute(fid, vid, 'units', len=attr_len)
      if( rcc.ne.nf90_noerr ) then
         write(*, '(a)') ' WARN::load_timevalues:'//&
              'attribute ---units--- not available for time variable!'
      else if( attr_len.gt.256 ) then
         write(*, '(a,(2a,i4,1x))') ' WARN::load_timevalues:'//&
              'will not read units attribute since length exceeds internal buffer size: ',&
              'limit=',256,'got=',attr_len
      else
         rcc = nf90_get_att(fid, vid, 'units', values=time_unit)
         if( rcc.ne.nf90_noerr ) then
            write(*, '(a)') ' WARN::load_timevalues:'//&
                 'error while reading ---units--- attribute.'
            time_unit = ''
         endif
      endif

      !-- successfully loaded
      load_timevalues = .true.

    end function load_timevalues

end subroutine ncrd_retrieval_config


!***********************************************************
!     ncrd_retrieval_prior
!
!> @brief Initialiases and sets all required prior information which is necessary to
!         run the observational operator.
!
!> @details TBD
!
!> @param[in]  fname        file name of NetCDF files providing prior information
!> @param[out] succeed      flag indicating whether prior setup was successfully read
!
subroutine ncrd_retrieval_prior(fname, succeed)
  use netcdf
  use mo_prior
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname
  logical, intent(out) :: succeed
  ! local decls
  integer :: fid,rc
  integer :: rd_npts, rd_nparam_s1
  real(kind=8) :: state_buf(nsc,npts)    ! nsc x npts
  real(kind=8) :: stateunc_buf(nsc,npts) ! nsc x npts
  real(kind=8) :: param_s1_buf(nparam_s1), paramunc_s1_buf(nparam_s1) ! nparam_s1
  integer :: xshape(1)

  !-- initialise output flag
  succeed = .false.

  !-- open file handle
  rc = nf90_open(trim(fname), nf90_nowrite, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config:opening file failed',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- consistent dimension of parameter for S1
  rd_nparam_s1 = get_dimsiz('nparam_s1')
  if( rd_nparam_s1.lt.0 ) then
     return
  else if( rd_nparam_s1.ne.nparam_s1 ) then
     return
  endif
  !-- read parameter for S1
  if( .not. read_param_s1() ) then
     return
  else
     
  endif

  !-- consistent dimension of statevector
  rd_npts = get_dimsiz('npoints')
  if( rd_npts.lt.0 ) then
     return
  else if( rd_npts.ne.npts ) then
     return
  else
  endif
  !-- read state vector
  if( .not. read_states() ) then
     return
  else

  endif

  !-- allocate 1D control vector and uncertainties
  allocate( xpr(nparam_s1+nsc*npts), sxpr(nparam_s1+nsc*npts) )

  !-- map parameter for S1
  xpr(1:nparam_s1)  = param_s1_buf
  sxpr(1:nparam_s1) = paramunc_s1_buf
  !-- map state vector
  xshape(1) = nsc*npts
  xpr(nparam_s1+1:)  = reshape( state_buf, xshape )
  sxpr(nparam_s1+1:) = reshape( stateunc_buf, xshape )


  !-- close handle
  rc = nf90_close(fid)

  !-- eventually indicate successful termination
  succeed = .true.


  contains

    integer function get_dimsiz(dim_name)
      implicit none
      character(len=*), intent(in) :: dim_name
      integer dimid,dimsiz,rcc

      get_dimsiz = -1
      rcc = nf90_inq_dimid(fid, dim_name, dimid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_prior::get_dimsiz:'//&
              'npoints dimension not accessible, ',&
              trim(nf90_strerror(rcc))
         return
      else
         rcc = nf90_inquire_dimension(fid, dimid, len=dimsiz)
         if( rcc.ne.nf90_noerr ) then
            write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config::get_dimsiz:'//&
                 'npoints dimension extent undeterminable',&
                 trim(nf90_strerror(rcc))
            return
         else
            get_dimsiz = dimsiz
         endif
      endif
    end function get_dimsiz

    logical function read_param_s1()
      implicit none
      integer :: rcc,vid

      !-- initialise succeed flag
      read_param_s1 = .false.

      !-- parameter:
      rcc = nf90_inq_varid(fid, 'param_s1', vid)
      rcc = nf90_get_var(fid, vid, param_s1_buf)
      
      !-- parameter uncertainty
      rcc = nf90_inq_varid(fid, 'param_unc_s1', vid)
      rcc = nf90_get_var(fid, vid, paramunc_s1_buf)

      !-- eventually set succeed flag
      read_param_s1 = .true.
    end function read_param_s1

    logical function read_states()
      implicit none
      integer :: rcc,vid

      read_states = .false.

      !-- LAI
      rcc = nf90_inq_varid(fid, 'lai', vid)
      rcc = nf90_get_var(fid, vid, state_buf(1,:))
      !-- Canopy Height
      rcc = nf90_inq_varid(fid, 'canht', vid)
      rcc = nf90_get_var(fid, vid, state_buf(2,:))
      !-- Soil moisture
      rcc = nf90_inq_varid(fid, 'sm', vid)
      rcc = nf90_get_var(fid, vid, state_buf(3,:))

      !-- LAI uncertainty
      rcc = nf90_inq_varid(fid, 'lai_unc', vid)
      rcc = nf90_get_var(fid, vid, stateunc_buf(1,:))
      !-- Canopy Height uncertainty
      rcc = nf90_inq_varid(fid, 'canht_unc', vid)
      rcc = nf90_get_var(fid, vid, stateunc_buf(2,:))
      !-- Soil moisture uncertainty
      rcc = nf90_inq_varid(fid, 'sm_unc', vid)
      rcc = nf90_get_var(fid, vid, stateunc_buf(3,:))

      read_states = .true.

    end function read_states
end subroutine ncrd_retrieval_prior


!***********************************************************
!     ncrd_retrieval_model
!
!> @brief reads (and sets) required environment to run the dynamic model
!
!> @details TBD
!
!> @param[in]  fname        file name of NetCDF files providing prior information
!> @param[out] succeed      flag indicating whether prior setup was successfully read
!
subroutine ncrd_retrieval_model(fname, succeed)
  use netcdf
  use mo_model
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname
  logical, intent(out) :: succeed
  ! local decls
  integer :: fid,rc
  integer :: rd_npts
  real(kind=8), allocatable :: a(:,:) !lai,canht,sm
  real(kind=8), allocatable :: b(:,:) !lai,canht,sm
  real(kind=8), allocatable :: munc(:,:) !
  integer :: i,j,k

  !-- initialise output flag
  succeed = .false.

  !-- open file handle
  rc = nf90_open(trim(fname), nf90_nowrite, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config:opening file failed',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- consistent dimension of model (in NetCDF file)
  rd_npts = get_dimsiz('npoints')
  if( rd_npts.lt.0 ) then
     return
  else if( rd_npts.ne.npts ) then
     write(*,'(a,2(a,i4,1x))') ' ERROR::ncrd_retrieval_model:'//&
          'unexpected dimension','expected=',npts,'got=',rd_npts
     return
  else
     allocate( a(npts,nsc), b(npts,nsc), munc(npts,nsc) )
  endif

  !-- read state vector
  if( .not. read_data() ) then
     return
  else
     a_model(1:nstates:nsc) = a(:,1)
     do i=1,npts
        do j=1,nsc
           k = (i-1)*nsc+j
           a_model(k)      = a(i,j)
           offset_model(k) = b(i,j)
           unc_model(k)    = munc(i,j)
        enddo
     enddo
  endif


  !-- dispose memory
  if( allocated(a) )    deallocate(a)
  if( allocated(b) )    deallocate(b)
  if( allocated(munc) ) deallocate(munc)

  !-- close handle
  rc = nf90_close(fid)

  !-- eventually indicate successful termination
  succeed = .true.

  contains

    integer function get_dimsiz(dim_name)
      implicit none
      character(len=*), intent(in) :: dim_name
      integer dimid,dimsiz,rcc

      get_dimsiz = -1
      rcc = nf90_inq_dimid(fid, dim_name, dimid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_prior::get_dimsiz:'//&
              'npoints dimension not accessible, ',&
              trim(nf90_strerror(rcc))
         return
      else
         rcc = nf90_inquire_dimension(fid, dimid, len=dimsiz)
         if( rcc.ne.nf90_noerr ) then
            write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_retrieval_config::get_dimsiz:'//&
                 'npoints dimension extent undeterminable',&
                 trim(nf90_strerror(rcc))
            return
         else
            get_dimsiz = dimsiz
         endif
      endif
    end function get_dimsiz

    logical function read_data()
      implicit none
      integer :: rcc,vid

      read_data = .false.

      !--     a - v a l u e s
      !-- LAI
      rcc = nf90_inq_varid(fid, 'a_lai', vid)
      rcc = nf90_get_var(fid, vid, a(:,1) )
      !-- Canopy Height
      rcc = nf90_inq_varid(fid, 'a_canht', vid)
      rcc = nf90_get_var(fid, vid, a(:,2))
      !-- Soil moisture
      rcc = nf90_inq_varid(fid, 'a_sm', vid)
      rcc = nf90_get_var(fid, vid, a(:,3))

      !--     b - v a l u e s
      !-- LAI
      rcc = nf90_inq_varid(fid, 'offset_lai', vid)
      rcc = nf90_get_var(fid, vid, b(:,1) )
      !-- Canopy Height
      rcc = nf90_inq_varid(fid, 'offset_canht', vid)
      rcc = nf90_get_var(fid, vid, b(:,2))
      !-- Soil moisture
      rcc = nf90_inq_varid(fid, 'offset_sm', vid)
      rcc = nf90_get_var(fid, vid, b(:,3))

      !--      m u n c - v a l u e s
      rcc = nf90_inq_varid(fid, 'munc_lai', vid)
      rcc = nf90_get_var(fid, vid, munc(:,1) )
      !-- Canopy Height
      rcc = nf90_inq_varid(fid, 'munc_canht', vid)
      rcc = nf90_get_var(fid, vid, munc(:,2))
      !-- Soil moisture
      rcc = nf90_inq_varid(fid, 'munc_sm', vid)
      rcc = nf90_get_var(fid, vid, munc(:,3))

      read_data = .true.

    end function read_data

end subroutine ncrd_retrieval_model


!***********************************************************
!     ncwrt_simulation_vector
!
!> @brief Writes Sentinel Simulator generated (1D) simulation vector for
!         specific simulation typ (S1 or S2) to NetCDF file.
!
!> @details TBD
!
!> @param[in]  fname        file name of NetCDF files providing prior information
!> @param[in]  simtyp       simulation typ, either 's1' or 's2'
!> @param[in]  m            1D-size of simulation vector (for current simulation typ)
!> @param[in]  y            1D simulation vector (for current simulation typ)
!> @param[out] succeed      flag indicating whether prior setup was successfully read
!
subroutine ncwrt_simulation_vector(fname, simtyp, m, y, yunc, succeed)
  use netcdf
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname
  character(len=2), intent(in) :: simtyp !'s1','s2'
  integer, intent(in) :: m
  real(kind=8), intent(in) :: y(m), yunc(m)
  logical, intent(out) :: succeed
  ! local decls
  integer :: rc, fid, data_vid, dataunc_vid
  integer :: time_vid
  integer :: dim_ids(2)
  integer :: ndim1, ndim2
  character(len=16) :: dataset_name, datasetunc_name, dim1, dim2
  real(kind=8), allocatable :: timepts_out(:)
  integer :: i,isim

  !-- initialise return flag
  succeed = .false.

  select case(simtyp)
     case('s1')
        ndim1 = get_npol()
        ndim2 = npts_s1
        dataset_name = 'backscatter'
        datasetunc_name = 'backscatter_unc'
        dim1 = 'npol'
     case('s2')
        ndim1 = get_nbs2()
        ndim2 = npts_s2
        dataset_name = 'brf'
        datasetunc_name = 'brf_unc'
        dim1 = 'nbands'
     case default
        write(*,'()') ' FATAL::ncwrt_simulation_vector:unexpected simtyp +++'//&
             trim(simtyp)//'+++ cannot continue.'
        return
  end select
  dim2 = 'npoints'

  !-- dimension consistency
  if( m.ne.ndim1*ndim2 ) then
     write(*,'(a,2(a,i3,1x))') ' FATAL::ncwrt_simulation_vector:'//&
          'inconsistent dimensions ','m=',m,'ndim1*ndim2=',ndim1*ndim2
     return
  endif

  !-- open file handle
  rc = nf90_create(trim(fname), NF90_NETCDF4, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector: opening file '//&
          '***'//trim(fname)//'*** for writing failed: ',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- get time-values
  allocate( timepts_out(ndim2) )
  isim = 1
  do i=1,npts
     if( simtyp.eq.'s1' .and. idx_is_s1(i) ) then
        timepts_out(isim) = timepts(i)
        isim = isim+1
     endif
     if( simtyp.eq.'s2' .and. idx_is_s2(i) ) then
        timepts_out(isim) = timepts(i)
        isim = isim+1
     endif
  enddo


  !-- create dimensions
  rc = nf90_def_dim(fid, dim1, ndim1, dim_ids(1))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector: creation of dimension '//&
          '---'//trim(dim1)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if
  rc = nf90_def_dim(fid, dim2, ndim2, dim_ids(2))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector: creation of dimension '//&
          '---'//trim(dim2)//'--- failed.', trim(nf90_strerror(rc))
     return
  end if

  !-- define time variable
  rc = nf90_def_var(fid, 'time', NF90_REAL8, dim_ids(2), time_vid, deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector:'//&
          'creation of dataset variable '//&
          '---'//trim('time')//'--- failed.', trim(nf90_strerror(rc))
     return
  endif

  !-- define dataset variable
  rc = nf90_def_var(fid, dataset_name, NF90_REAL8, dim_ids, data_vid, deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name)//'--- failed.', trim(nf90_strerror(rc))
     return
  endif

  !-- define dataset uncertainty variable
  rc = nf90_def_var(fid, datasetunc_name, NF90_REAL8, dim_ids, dataunc_vid, deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector:'//&
          'creation of dataset uncertainty variable '//&
          '---'//trim(datasetunc_name)//'--- failed.', trim(nf90_strerror(rc))
     return
  endif

  !-------------------
  !     e n d   d e f i n e   m o d e
  !
  rc = nf90_enddef(fid)

  !-- write time
  rc = nf90_put_var(fid, time_vid, timepts_out)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector:'//&
          'error when writing data to file! ', trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, time_vid, 'standard_name', 'time')
  rc = nf90_put_att(fid, time_vid, 'long_name', 'time')
  if( time_unit.ne.'' ) then
     rc = nf90_put_att(fid, time_vid, 'units', time_unit)
  endif

  !-- write data
  rc = nf90_put_var(fid, data_vid, reshape(y, (/ ndim1, ndim2 /)))
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector:'//&
          'error when writing data to file! ', trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, data_vid, 'missing_value', sim_fill_value)

  !-- write data uncertainty
  rc = nf90_put_var(fid, dataunc_vid, reshape(yunc, (/ ndim1, ndim2 /)))
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_simulation_vector:'//&
          'error when writing data to file! ', trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, dataunc_vid, 'missing_value', sim_fill_value)

  !-- close file
  rc = nf90_close(fid)

  !-- dispose memory
  if( allocated(timepts_out) ) deallocate(timepts_out)

  !-- 
  succeed = .true.

end subroutine ncwrt_simulation_vector


!***********************************************************
!     ncwrt_target_vector
!
!> @brief Writes Sentinel Simulator generated (1D) simulation vector for
!         specific simulation typ (S1 or S2) to NetCDF file.
!
!> @details TBD
!
!> @param[in]  fname        file name of NetCDF files providing prior information
!> @param[in]  dsname       dataset identifier
!> @param[in]  m            1D-size of simulation vector (for current simulation typ)
!> @param[in]  y            1D simulation vector (for current simulation typ)
!> @param[out] succeed      flag indicating whether prior setup was successfully read
!
subroutine ncwrt_target_vector(fname, dsname, m, y, yunc, succeed)
  use netcdf
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: dsname
  integer, intent(in) :: m
  real(kind=8), intent(in) :: y(m), yunc(m)
  logical, intent(out) :: succeed
  ! local decls
  integer :: rc, fid, data_vid, dataunc_vid
  integer :: time_vid
  integer :: dim_ids(1)
  character(len=16) :: dataset_name, datasetunc_name, dim1, dim2
  real(kind=8), allocatable :: timepts_out(:)

  !-- initialise return flag
  succeed = .false.

  if( m.ne.get_npts() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::ncwrt_target_vectorr:'//&
          'inconsistent length of simulation vector.',&
          'expected=',get_npts(),'got=',m
  endif
  

  dataset_name = dsname
  datasetunc_name = trim(dsname)//'_unc'
  dim2 = 'npoints'


  !-- open file handle
  rc = nf90_create(trim(fname), NF90_NETCDF4, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector: opening file '//&
          '***'//trim(fname)//'*** for writing failed: ',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- get time-values
  allocate( timepts_out(m) )
  timepts_out = timepts


  !-- create dimensions
  rc = nf90_def_dim(fid, 'npoints', m, dim_ids(1))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector: creation of dimension '//&
          '---'//trim(dim1)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if

  !-- define time variable
  rc = nf90_def_var(fid, 'time', NF90_REAL8, dim_ids(1), time_vid, deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector:'//&
          'creation of dataset variable '//&
          '---'//trim('time')//'--- failed.', trim(nf90_strerror(rc))
     return
  endif

  !-- define dataset variable
  rc = nf90_def_var(fid, dataset_name, NF90_REAL8, dim_ids, data_vid, deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name)//'--- failed.', trim(nf90_strerror(rc))
     return
  endif

  !-- define dataset uncertainty variable
  rc = nf90_def_var(fid, datasetunc_name, NF90_REAL8, dim_ids, dataunc_vid, deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector:'//&
          'creation of dataset uncertainty variable '//&
          '---'//trim(datasetunc_name)//'--- failed.', trim(nf90_strerror(rc))
     return
  endif

  !-------------------
  !     e n d   d e f i n e   m o d e
  !
  rc = nf90_enddef(fid)

  !-- write time
  rc = nf90_put_var(fid, time_vid, timepts_out)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector:'//&
          'error when writing data to file! ', trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, time_vid, 'standard_name', 'time')
  rc = nf90_put_att(fid, time_vid, 'long_name', 'time')
  if( time_unit.ne.'' ) then
     rc = nf90_put_att(fid, time_vid, 'units', time_unit)
  endif

  !-- write data
  rc = nf90_put_var(fid, data_vid, y)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector:'//&
          'error when writing data to file! ', trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, data_vid, 'missing_value', sim_fill_value)

  !-- write data uncertainty
  rc = nf90_put_var(fid, dataunc_vid, yunc)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_target_vector:'//&
          'error when writing data to file! ', trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, dataunc_vid, 'missing_value', sim_fill_value)

  !-- close file
  rc = nf90_close(fid)

  !-- dispose memory
  if( allocated(timepts_out) ) deallocate(timepts_out)

  !-- 
  succeed = .true.

end subroutine ncwrt_target_vector


!***********************************************************
!     ncrd_obs_vector
!
!> @brief Reades observation vector and associated unceratainties from NetCDF file.
!
!> @details TBD
!
!> @param[in]  fname        file name of NetCDF files providing prior information
!> @param[in]  obstyp       observation typ, either 's1' or 's2'
!> @param[in]  m            1D-size of observation vector
!> @param[in]  y            1D observation vector (for current obstyp)
!> @param[out] succeed      flag indicating whether prior setup was successfully read
!
subroutine ncrd_obs_vector(fname, obstyp, m, obsvec, obsuncvec, succeed)
  use netcdf
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname
  character(len=2), intent(in) :: obstyp !'s1','s2'
  integer, intent(in) :: m
  real(kind=8), intent(out) :: obsvec(m), obsuncvec(m)
  logical, intent(out) :: succeed
  ! local decls
  integer :: rc, fid, data_vid, dataunc_vid
  integer :: cnt(1)
  integer :: ndim1, ndim2
  character(len=16) :: dataset_name, datasetunc_name, dim1, dim2
  real(kind=8), allocatable :: io_buffer(:,:)
  real(kind=8) :: data_fv, dataunc_fv
  !-- initialise return flag
  succeed = .false.

  select case(obstyp)
     case('s1')
        ndim1 = get_npol()
        ndim2 = npts_s1
        dataset_name = 'backscatter'
        datasetunc_name = 'backscatter_unc'
        dim1 = 'npol'
     case('s2')
        ndim1 = get_nbs2()
        ndim2 = npts_s2
        dataset_name = 'brf'
        datasetunc_name = 'brf_unc'
        dim1 = 'nbands'
     case default
        write(*,'()') ' FATAL::ncrd_obs_vector:unexpected obstyp +++'//&
             trim(obstyp)//'+++ cannot continue.'
        return
  end select
  dim2 = 'npoints'

  !-- dimension consistency
  if( m.ne.ndim1*ndim2 ) then
     write(*,'(a,2(a,i3,1x))') ' FATAL::ncrd_obs_vector:'//&
          'inconsistent dimensions ','m=',m,'ndim1*ndim2=',ndim1*ndim2
     return
  endif

  !-- open file handle
  rc = nf90_open(trim(fname), nf90_nowrite, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_obs_vector:'//&
          'opening file ***'//trim(fname)//'*** failed',&
          trim(nf90_strerror(rc))
     return
  end if

  if( .not. check_dim_sizes() ) then
     write(*,'(a)') ' ERROR::ncrd_obs_vector:dimension-check failed, observations are not loaded.'
     return
  endif

  !-- determine handles: data set
  rc = nf90_inq_varid(fid, dataset_name, data_vid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_obs_vector:'//&
          'dataset ---'//trim(dataset_name)//'--- not found in file: ',&
          trim(nf90_strerror(rc))
     return
  end if

  ! inquire fill-value
  rc = nf90_get_att(fid, data_vid, 'missing_value', values=data_fv)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::',trim(nf90_strerror(rc))
     return
  else
     write(*,'(a,f15.8)') ' INFO::ncrd_obs_vector:'//&
          'dataset ---'//trim(dataset_name)//'--- detected fill-value:',data_fv
  end if


  !-- determine handles: data set uncertainty
  rc = nf90_inq_varid(fid, datasetunc_name, dataunc_vid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_obs_vector:'//&
          'dataset uncertainty ---'//trim(datasetunc_name)//'--- not found in file: ',&
          trim(nf90_strerror(rc))
     return
  end if

  ! inquire fill-value
  rc = nf90_get_att(fid, dataunc_vid, 'missing_value', values=dataunc_fv)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::',trim(nf90_strerror(rc))
     return
  else
     write(*,'(a,f15.8)') ' INFO::ncrd_obs_vector:'//&
          'dataset ---'//trim(datasetunc_name)//'--- detected fill-value:',dataunc_fv
  end if

  !-- allocate I/O buffer
  allocate( io_buffer(ndim1,ndim2) )

  !-- read data set
  rc = nf90_get_var( fid, data_vid, io_buffer )
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_obs_vector:'//&
          'error while reading dataset ---'//trim(dataset_name)//'--- : ',&
          trim(nf90_strerror(rc))
     return
  else
     obsvec = reshape(io_buffer, (/m/))
  end if

  !-- remap obs. fill-values to assimilation fill-values
  cnt = count(obsvec.eq.data_fv)
  if( cnt(1).gt.0 ) then
     write(*, '(a,i8)' ) ' INFO::ncrd_obs_vector:'//&
          obstyp//' detected data fill-values:',cnt(1)
     where( obsvec.eq.data_fv )
        obsvec = sim_fill_value
     end where
  else
     write(*, '(a)' ) ' INFO::ncrd_obs_vector:'//&
          obstyp//' detected no fill-values'
  endif
  !-- 
  where( obsvec.eq.data_fv )
     obsvec = sim_fill_value
  end where

  !-- read data set uncertainty
  io_buffer = -1._8
  rc = nf90_get_var( fid, dataunc_vid, io_buffer )
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncrd_obs_vector:'//&
          'error while reading dataset ---'//trim(dataset_name)//'--- : ',&
          trim(nf90_strerror(rc))
     return
  else
     obsuncvec = reshape(io_buffer, (/m/))
  end if

  !-- remap unc. fill-values to assimilation fill-values
  cnt = count(obsuncvec.eq.dataunc_fv)
  !-- 
  if( cnt(1).gt.0 ) then
     write(*, '(a,i8)' ) ' INFO::ncrd_obs_vector:'//&
          obstyp//' detected unc. fill-values:',cnt(1)
     where( obsuncvec.eq.dataunc_fv )
        obsuncvec = sim_fill_value
     end where
  else
     write(*, '(a)' ) ' INFO::ncrd_obs_vector:'//&
          obstyp//' detected no uncertainthy fill-values'
  endif

  !-- merge fill-value ranges (obs. and obs unc.)
  where( obsvec.eq.sim_fill_value.or.obsuncvec.eq.sim_fill_value )
     obsvec    = sim_fill_value
     obsuncvec = sim_fill_value
  end where

  !-- dispose memory
  deallocate( io_buffer )

  !-- successfully read data
  succeed = .true.

  !-- close handle
  rc = nf90_close(fid)

  contains
    logical function check_dim_sizes()
      implicit none
      integer :: rcc,dimid,dimsiz

      check_dim_sizes = .false.

      !-- check first dimension
      rcc = nf90_inq_dimid(fid, dim1, dimid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::check_dim_sizes: dimension='//trim(dim1)//': ',&
              trim(nf90_strerror(rcc))
         return
      else
         rcc = nf90_inquire_dimension(fid, dimid, len=dimsiz)
         if( rcc.ne.nf90_noerr ) then
            write(*,'(a,a)') ' ###NETCDF-ERROR::check_dim_sizes: dimension extent',&
                 trim(nf90_strerror(rcc))
            return
         else if( dimsiz.ne.ndim1 ) then
            write(*,'(a,2(a,i3,1x))') ' ###NETCDF-ERROR::check_dim_sizes: dimension='//trim(dim1)//&
                 ' with unexpected size, ','expected=',ndim1,'got=',dimsiz
            return
         endif
      endif

      !-- check second dimension
      rcc = nf90_inq_dimid(fid, dim2, dimid)
      if( rcc.ne.nf90_noerr) then
         write(*,'(a,a)') ' ###NETCDF-ERROR::check_dim_sizes: dimension='//trim(dim2)//': ',&
              trim(nf90_strerror(rcc))
         return
      else
         rcc = nf90_inquire_dimension(fid, dimid, len=dimsiz)
         if( rcc.ne.nf90_noerr ) then
            write(*,'(a,a)') ' ###NETCDF-ERROR::check_dim_sizes: dimension extent',&
                 trim(nf90_strerror(rcc))
            return
         else if( dimsiz.ne.ndim2 ) then
            write(*,'(a,2(a,i3,1x))') ' ###NETCDF-ERROR::check_dim_sizes: dimension='//trim(dim2)//&
                 ' with unexpected size, ','expected=',ndim2,'got=',dimsiz
            return
         endif
      endif

      !-- dimension check succeeded
      check_dim_sizes = .true.
    end function check_dim_sizes
end subroutine ncrd_obs_vector


!***********************************************************
!     ncwrt_state_and_sim
!
!> @brief Postprocessing routine to consistently write
!>        control vector and simulation vector to NetCDF file.
!
subroutine ncwrt_state_and_sim(fname, n, x, sx, m, y, succeed)
  use netcdf
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in) :: fname   !< file name of NetCDF to be generated
  integer, intent(in)          :: n       !< length of control vector
  integer, intent(in)          :: m       !< length of simulation vector
  real(kind=8), intent(in)     :: x(n)    !< control vector
  real(kind=8), intent(in)     :: sx(n)   !< prior uncertainty of control vector
  real(kind=8), intent(in)     :: y(m)    !< simulation vector
  logical, intent(out)         :: succeed !< flag indicating whether I/O succeeded
  ! local decls
  integer, parameter :: nvar = nsc+nsc+1+1+1+1 !time statecomponents,their uncertainties,hhvv,brf,timept_simtyp
  real(kind=8), parameter :: state_fill_value = -99999.
  real(kind=8), allocatable :: statev(:,:), stateuncv(:,:), brf(:,:), hhvv(:,:)
  real(kind=8), allocatable :: brf_buf(:,:), hhvv_buf(:,:)
  integer :: np, m_s1, m_s2
  integer :: i,i_s1,i_s2,ivar
  integer :: rc, fid
  integer :: dataset_ids(nvar)
  character(len=16) :: dataset_name(nvar)
  character(len=16) :: dim_name
  integer :: dim_ids(3)  !npol,nb_s2,npts
  integer :: npol, nb_s2

  !-- initialise return flag
  succeed = .false.

  !-- configuration sizes
  npol = get_npol()
  nb_s2 = get_nbs2()
  np = get_np()
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()

  !-- remap 1D state vector
  allocate( statev(nsc,npts), stateuncv(nsc,npts), brf(nb_s2,npts), hhvv(npol,npts) )
  statev    = reshape( x(np+1:),  (/nsc,npts/) )
  stateuncv = reshape( sx(np+1:), (/nsc,npts/) )
  allocate(brf_buf(nb_s2,npts),hhvv_buf(npol,npts))
  hhvv_buf = reshape( y(1:m_s1), (/npol,npts_s1/) )
  brf_buf  = reshape( y(m_s1+1:m_s1+m_s2), (/nb_s2,npts_s2/) )
  hhvv = state_fill_value
  brf  = state_fill_value
  i_s1 = 1
  i_s2 = 1
  do i=1,npts
     if( idx_is_s1(i) ) then
        hhvv(:,i) = hhvv_buf(:,i_s1)
        i_s1 = i_s1 + 1
     endif
     if( idx_is_s2(i) ) then
        brf(:,i)  = brf_buf(:,i_s2)
        i_s2 = i_s2 + 1
     endif
  enddo
  deallocate(hhvv_buf, brf_buf)

  !-- open file handle
  rc = nf90_create(trim(fname), NF90_NETCDF4, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim: opening file '//&
          '***'//trim(fname)//'*** for writing failed: ',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- create dimensions
  dim_name = 'npol'
  rc = nf90_def_dim(fid, dim_name, npol, dim_ids(1))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim: creation of dimension '//&
          '---'//trim(dim_name)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if

  dim_name = 'nb_s2'
  rc = nf90_def_dim(fid, dim_name, nb_s2, dim_ids(2))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim: creation of dimension '//&
          '---'//trim(dim_name)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if
  dim_name = 'npoints'
  rc = nf90_def_dim(fid, dim_name, npts, dim_ids(3))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim: creation of dimension '//&
          '---'//trim(dim_name)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if

  !-- define dataset variables
  ivar = 1

  dataset_name(ivar) = 'time'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'time_typ'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_INT,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'lai'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'canht'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'sm'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'lai_unc'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'canht_unc'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'sm_unc'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       dim_ids(3), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'backscatter'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       (/dim_ids(1),dim_ids(3)/), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  dataset_name(ivar) = 'brf'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
       (/dim_ids(2),dim_ids(3)/), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar+1

  !--     e n d   d e f i n e   m o d e
  rc = nf90_enddef(fid)

  !-- start writing data...
  ivar = 1

  !-- time
  rc = nf90_put_var(fid, dataset_ids(ivar), timepts)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
          'error when writing data for var ---'//trim(dataset_name(ivar))//'---   ',&
          trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'time')
  rc = nf90_put_att(fid, dataset_ids(ivar), 'long_name', 'time')
  if( time_unit.ne.'' ) then
     rc = nf90_put_att(fid, dataset_ids(ivar), 'units', time_unit)
  endif
  ivar = ivar + 1

  !-- time_typ
  rc = nf90_put_var(fid, dataset_ids(ivar), timept_simtyp)
  if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
             'error when writing data for var ---'//trim(dataset_name(ivar))//'---   ',&
             trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  !-- states
  do i=1,nsc
     call set_attributes(ivar)
     rc = nf90_put_var(fid, dataset_ids(ivar), statev(i,:))
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
             'error when writing data for var ---'//trim(dataset_name(ivar))//'---   ',&
             trim(nf90_strerror(rc))
        return
     endif
     ivar = ivar+1
  enddo

  !-- prior state uncertainty
  do i=1,nsc
     call set_attributes(ivar)
     rc = nf90_put_var(fid, dataset_ids(ivar), stateuncv(i,:))
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
             'error when writing data for var ---'//trim(dataset_name(ivar))//'---   ',&
             trim(nf90_strerror(rc))
        return
     endif
     ivar = ivar+1
  enddo

  !-- S1 simulations
  call set_attributes(ivar)
  rc = nf90_put_var(fid, dataset_ids(ivar), hhvv(:,:))
  if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
             'error when writing data for var ---'//trim(dataset_name(ivar))//'---   ',&
             trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  !-- S2 simulations
  call set_attributes(ivar)
  rc = nf90_put_var(fid, dataset_ids(ivar), brf(:,:))
  if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_state_and_sim:'//&
             'error when writing data for var ---'//trim(dataset_name(ivar))//'---   ',&
             trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  ! !-- add global attributes
  ! call set_global_attributes()

  !-- close file
  rc = nf90_close(fid)

  !-- successfully written the file
  succeed = .true.

  contains
    subroutine set_attributes(ivar)
      implicit none
      integer, intent(in) :: ivar
      ! local decls
      integer :: rcc

      rcc = nf90_put_att(fid, dataset_ids(ivar), 'missing_value', state_fill_value)

      select case( dataset_name(ivar) )
         case('lai')
            rcc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'LAI')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'long_name', 'Leaf area index')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'units','m2/m2')
         case('canht')
            rcc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'canht')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'long_name', 'Canopy height')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'units','m')
         case('sm')
            rcc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'SM')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'long_name', 'Soil moisture')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'units','m3/m3')
         case('lai_unc')
            rcc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'LAI_unc')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'long_name', '1-sigma uncertainty of Leaf area index')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'units','m2/m2')
         case('canht_unc')
            rcc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'canht_unc')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'long_name', '1-sigma uncertainty of Canopy height')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'units','m')
         case('sm_unc')
            rcc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'SM_unc')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'long_name', '1-sigma uncertainty of Soil moisture')
            rcc = nf90_put_att(fid, dataset_ids(ivar),&
                 'units','m3/m3')
         case default
            continue
      end select
    end subroutine set_attributes

    subroutine set_global_attributes()
      use mo_timing
      implicit none
      integer :: rcc
      rcc = nf90_put_att(fid, nf90_global,&
           'institution', 'The Inversion Lab, Hamburg, Germany')
      rcc = nf90_put_att(fid, nf90_global,&
           'date_created', time_now())
    end subroutine set_global_attributes

end subroutine ncwrt_state_and_sim


!***********************************************************
!     ncwrt_controlvector
!
!> @brief Postprocessing routine to consistently write
!>        control vector to NetCDF file.
!
logical function ncwrt_controlvector(fname, n, x, sx, optiflag) result(succeed)
  use netcdf
  use mo_sensimul
  implicit none
  ! arguments
  character(len=*), intent(in)  :: fname
  integer, intent(in)           :: n
  real(kind=8), intent(in)      :: x(n)
  real(kind=8), intent(in)      :: sx(n)
  logical, intent(in), optional :: optiflag
  ! local decls
  integer, parameter :: nvar = 1+1+2*(1+nsc) !time,state-typ,
                                             !parameter,parameter uncertainty, statecomponents,their uncertainties,
  real(kind=8), parameter :: state_fill_value = -99999.
  real(kind=8), allocatable :: ctlvec_params(:), ctlvec_params_unc(:)
  real(kind=8), allocatable :: statev(:,:), stateuncv(:,:)
  integer :: np
  integer :: rc, fid, dim_id(2)
  integer :: dataset_ids(nvar) 
  character(len=16) :: dataset_name(nvar)
  character(len=16) :: dim_name
  integer :: i,ivar

  !-- initialise return flag
  succeed = .false.

  !-- remap 1D state vector
  np = get_np()
  if( np.ne.1 ) then
     write(*, '(a)') ' FATAL::ncwrt_controlvector:actually only nparam=1 is supported!'
     stop
  endif

  allocate( ctlvec_params(np), ctlvec_params_unc(np) )
  allocate( statev(nsc,npts), stateuncv(nsc,npts) )

  !--
  ctlvec_params(1:np)     = x(1:np)
  ctlvec_params_unc(1:np) = sx(1:np)
  statev    = reshape( x(np+1:),  (/nsc,npts/) )
  stateuncv = reshape( sx(np+1:), (/nsc,npts/) )

  !-- open file handle
  rc = nf90_create(trim(fname), NF90_NETCDF4, fid)
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector: opening file '//&
          '***'//trim(fname)//'*** for writing failed: ',&
          trim(nf90_strerror(rc))
     return
  end if

  !-- create dimension(s)
  dim_name = 'nparam_s1'
  rc = nf90_def_dim(fid, dim_name, np, dim_id(1))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector: creation of dimension '//&
          '---'//trim(dim_name)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if

  dim_name = 'npoints'
  rc = nf90_def_dim(fid, dim_name, npts, dim_id(2))
  if( rc.ne.nf90_noerr) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector: creation of dimension '//&
          '---'//trim(dim_name)//'--- failed: ', trim(nf90_strerror(rc))
     return
  end if

  !-- define dataset variables
  ivar = 1

  !-- time
  dataset_name(ivar) = 'time'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8, dim_id(2), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  !-- time_typ
  dataset_name(ivar) = 'sim_typ'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_INT,&
       dim_id(2), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  !-- parameter
  dataset_name(ivar) = 'laicoeff'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8, dim_id(1), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  !-- parameter uncertainty
  dataset_name(ivar) = 'laicoeff_unc'
  rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8, dim_id(1), dataset_ids(ivar), deflate_level=4)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'creation of dataset variable '//&
          '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
     return
  endif
  ivar = ivar + 1

  !-- state vector
  do i=1,nsc
     select case(i)
     case(1)
        dataset_name(ivar) = 'lai'
     case(2)
        dataset_name(ivar) = 'canht'
     case(3)
        dataset_name(ivar) = 'sm'
     end select
     rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
          dim_id(2), dataset_ids(ivar), deflate_level=4)
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
             'creation of dataset variable '//&
             '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
        return
     endif
     ivar = ivar + 1
  enddo

  !-- state vector uncertainty
  do i=1,nsc
     select case(i)
     case(1)
        dataset_name(ivar) = 'lai_unc'
     case(2)
        dataset_name(ivar) = 'canht_unc'
     case(3)
        dataset_name(ivar) = 'sm_unc'
     end select
     rc = nf90_def_var(fid, dataset_name(ivar), NF90_REAL8,&
          dim_id(2), dataset_ids(ivar), deflate_level=4)
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
             'creation of dataset variable '//&
             '---'//trim(dataset_name(ivar))//'--- failed.', trim(nf90_strerror(rc))
        return
     endif
     ivar = ivar + 1
  enddo

  !-- e n d   d e f i n e   m o d e
  rc = nf90_enddef(fid)

  !-- start writing data
  ivar = 1

  !-- time variable
  rc = nf90_put_var(fid, dataset_ids(ivar), timepts)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'error when writing data to file ---'//trim(dataset_name(ivar))//'---   ',&
          trim(nf90_strerror(rc))
     return
  endif
  rc = nf90_put_att(fid, dataset_ids(ivar), 'standard_name', 'time')
  rc = nf90_put_att(fid, dataset_ids(ivar), 'long_name', 'time')
  if( time_unit.ne.'' ) then
     rc = nf90_put_att(fid, dataset_ids(ivar), 'units', time_unit)
  endif
  ivar = ivar + 1

  !-- time_type variable
  rc = nf90_put_var(fid, dataset_ids(ivar), timept_simtyp)
  if( rc.ne.nf90_noerr ) then
     write(*,'(a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'error when writing data to file ---'//trim(dataset_name(ivar))//'---   '//&
          '(###'//trim(nf90_strerror(rc))//'###)'
     return
  endif
  rc = nf90_put_att(fid, dataset_ids(ivar), 'long_name', 'state type')
  rc = nf90_put_att(fid, dataset_ids(ivar), 'comment', 'integer value  to be bit-interpreted')
  rc = nf90_put_att(fid, dataset_ids(ivar), 'nobits_set',  'time-point with other state')
  rc = nf90_put_att(fid, dataset_ids(ivar), 'bit0_is_set', 'time-point for S1 simulation')
  rc = nf90_put_att(fid, dataset_ids(ivar), 'bit1_is_set', 'time-point for S2 simulation')
  ivar = ivar + 1

  !-- parameter
  rc = nf90_put_var(fid, dataset_ids(ivar), ctlvec_params(1))
  if( rc.ne.nf90_noerr ) then
     write(*,'(a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'error when writing data to file ---'//trim(dataset_name(ivar))//'---   '//&
          '(###'//trim(nf90_strerror(rc))//'###)'
     return
  endif
  rc = nf90_put_att(fid, dataset_ids(ivar), 'long_name', &
       'LAI conversion coefficient (optical to microwave)')
  ivar = ivar + 1
  
  !-- parameter uncertainty
  rc = nf90_put_var(fid, dataset_ids(ivar), ctlvec_params_unc(1))
  if( rc.ne.nf90_noerr ) then
     write(*,'(a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
          'error when writing data to file ---'//trim(dataset_name(ivar))//'---   '//&
          '(###'//trim(nf90_strerror(rc))//'###)'
     return
  endif
  rc = nf90_put_att(fid, dataset_ids(ivar), 'long_name', &
       '1-sigma uncertainty of LAI conversion coefficient (optical to microwave)')
  ivar = ivar + 1

  !-- state vector
  do i=1,nsc
     call set_attributes(ivar)
     rc = nf90_put_var(fid, dataset_ids(ivar), statev(i,:))
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
             'error when writing data to file for ---'//trim(dataset_name(ivar))//'--- ', trim(nf90_strerror(rc))
        return
     endif
     ivar = ivar + 1
  enddo

  !-- state vector uncertainty
  do i=1,nsc
     call set_attributes(ivar)
     rc = nf90_put_var(fid, dataset_ids(ivar), stateuncv(i,:))
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
             'error when writing data to file for ---'//trim(dataset_name(ivar))//'--- ', trim(nf90_strerror(rc))
        return
     endif
     ivar = ivar + 1
  enddo

  !-- add global attributes
  call set_global_attributes()

  !-- potentially add optimisation flag
  if( present(optiflag) ) then
     if( optiflag ) then
        rc = nf90_put_att(fid, nf90_global,&
             'optimisation_tag', 'SUCCESSFUL')
     else
        rc = nf90_put_att(fid, nf90_global,&
             'optimisation_tag', 'PROBLEMATIC')
     endif
     if( rc.ne.nf90_noerr ) then
        write(*,'(a,a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
             'error when writing optimisation flag to file ---'//trim(fname)//'--- ', trim(nf90_strerror(rc))
     endif
  endif

  !-- close file
  rc = nf90_close(fid)

  !-- 
  deallocate( ctlvec_params, ctlvec_params_unc )
  deallocate( statev, stateuncv )

  !-- successfully terminated I/O
  succeed = .true.



  contains
    subroutine set_attributes(ivar)
      implicit none
      ! arguments
      integer, intent(in) :: ivar
      ! local decls
      character(len=32) :: attname
      integer :: rcc

      attname = 'missing_value'
      rcc = nf90_put_att(fid, dataset_ids(ivar), attname, state_fill_value)
      if( rcc.ne.nf90_noerr ) then
         write(*,'(a)') ' ###NETCDF-ERROR::ncwrt_controlvector:'//&
             'error when setting attribute '//trim(attname)//&
             ' for ---'//trim(dataset_name(ivar))//'--- '//&
             '(###'//trim(nf90_strerror(rcc))//'###)'
        return
      endif

      select case( dataset_name(ivar) )
         case('lai')
            call set_cfvar_attributes(fid, dataset_ids(ivar), dataset_name(ivar),&
                 'LAI', 'Leaf area index', 'm2/m2')
         case('canht')
            call set_cfvar_attributes(fid, dataset_ids(ivar), dataset_name(ivar),&
                 'canht', 'Canopy height', 'm')
         case('sm')
            call set_cfvar_attributes(fid, dataset_ids(ivar), dataset_name(ivar),&
                 'SM', 'Soil moisture', 'm3/m3')
         case('lai_unc')
            call set_cfvar_attributes(fid, dataset_ids(ivar), dataset_name(ivar),&
                 'lai_unc', '1-sigma uncertainty of Leaf area index', 'm2/m2')
         case('canht_unc')
            call set_cfvar_attributes(fid, dataset_ids(ivar), dataset_name(ivar),&
                 'canht_unc', '1-sigma uncertainty of Canopy height', 'm')
         case('sm_unc')
            call set_cfvar_attributes(fid, dataset_ids(ivar), dataset_name(ivar),&
                 'sm_unc', '1-sigma uncertainty of Soil moisture', 'm3/m3')
         case default
            write(*, '(a)') ' FATAL::set_attributes:unexpected dataset ---'//&
                 trim(dataset_name(ivar))//'---'
            stop
      end select
    end subroutine set_attributes

    subroutine set_cfvar_attributes(fid, vid, varname, standard_name, long_name, units)
      implicit none
      integer, intent(in) :: fid, vid
      character(len=*), intent(in) :: varname, standard_name, long_name, units
      ! local declarations
      integer :: rcc
      character(len=32) :: attname

      attname = 'standard_name'
      rcc = nf90_put_att(fid, vid, attname, standard_name)
      if( rcc.ne.nf90_noerr ) then
         write(*, '(a)') ' NETCDF-ERROR::setting attribute *'//trim(attname)//&
              '* at variable ---'//trim(varname)//'--- did not succeed (###'//&
               trim(nf90_strerror(rcc))//'###)'
      endif

      attname = 'long_name'
      rcc = nf90_put_att(fid, vid, attname, long_name)
      if( rcc.ne.nf90_noerr ) then
         write(*, '(a)') ' NETCDF-ERROR::setting attribute *'//trim(attname)//&
              '* at variable ---'//trim(varname)//'--- did not succeed (###'//&
               trim(nf90_strerror(rcc))//'###)'
      endif

      attname = 'units'
      rcc = nf90_put_att(fid, vid, attname, units)
      if( rcc.ne.nf90_noerr ) then
         write(*, '(a)') ' NETCDF-ERROR::setting attribute *'//trim(attname)//&
              '* at variable ---'//trim(varname)//'--- did not succeed (###'//&
               trim(nf90_strerror(rcc))//'###)'
      endif
    end subroutine set_cfvar_attributes

    subroutine set_global_attributes()
      use mo_timing
      implicit none
      integer :: rcc
      rcc = nf90_put_att(fid, nf90_global,&
           'institution', 'The Inversion Lab, Hamburg, Germany')
      rcc = nf90_put_att(fid, nf90_global,&
           'date_created', time_now())
    end subroutine set_global_attributes

end function ncwrt_controlvector

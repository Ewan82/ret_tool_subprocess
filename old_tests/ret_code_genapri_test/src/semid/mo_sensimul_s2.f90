!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file mo_sensimul_s2.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: provision of variables used to run the optical model of the Sentinel Simulator (which are not regarded as "active" in the context of the retrieval tool.
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February, August 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mo_sensimul_s2
  implicit none

  private :: srf_file, read_s2_response

  ! wavelength range [400nm,2500nm] at 1nm resolution
  integer, parameter :: nw1nm = 2500 - 400 + 1 ! ==2101

  ! parameter NADIM
  real :: rpl = 0.1   ! radius of single leaf
  integer :: lad = 2  ! leaf-angle distribution
                      ! 1 - Planophile
                      ! 2 - Erectophile
                      ! 3 - Plagiophile
                      ! 4 - Extremophile
                      ! 5 - Uniform

  ! parameter PROSPECT
  real(kind=8) :: vai =  1.0   ! leaf structure
  real(kind=8) :: cab = 75.0   ! leaf chlorophyll
  real(kind=8) :: cw  =  0.01  ! leaf water equivalent thickness
  real(kind=8) :: cp  =  0.001 ! protein concentration
  real(kind=8) :: cc  =  0.001 ! cellulose and lignin

  ! parameter Price,
  ! rslN is the weight of the Nth spectral vector of the soil reflectance (of Price's EOF)
  ! IMPL-NOTE::rsl1 will be modulated by the Soil-Moisture state and thus become
  !            **active** in the sense of the retrieval system.
  real :: rsl1_default = 0.2
  real :: rsl2 =  0.1
  real :: rsl3 =  0.03726
  real :: rsl4 = -0.002426

  !-- weighting of soil moisture impact on rsl1, bound between (0,1)
  real :: sm_coeff = 0.5

  ! number of Sentinel2 wavebands:
  integer, parameter :: nb_s2 = 13  ! S2 bands B1-B12 (incl B8a)
  ! Sentinel2 wavebands:
  ! 1       443     Aerosole        20      60
  ! 2       490     Aerosole, Landnutzung, Vegetation       65      10
  ! 3       560     Landnutzung, Vegetation 35      10
  ! 4       665     Landnutzung, Vegetation 30      10
  ! 5       705     Landnutzung, Vegetation 15      20
  ! 6       740     Landnutzung, Vegetation 15      20
  ! 7       783     Landnutzung, Vegetation 20      20
  ! 8       842     Wasserdampf, Landnutzung, Vegetation    115     10
  ! 8a      865     Wasserdampf, Landnutzung, Vegetation    20      20
  ! 9       940     Wasserdampf     20      60
  ! 10      1375    Cirruswolken    20      60
  ! 11      1610    Landnutzung, Vegetation 90      20
  ! 12      2190    Aerosole, Landnutzung, Vegetation       180     20
  real(kind=8) :: wvl_s2(nb_S2) = (/&
       443._8, 490._8, 560._8, 665._8, 705._8,&
       740._8, 783._8, 842._8, 865._8, 940._8,&
       1375._8, 1610._8, 2190._8/)

  ! S2 responses ([400nm,2500nm] x 13 BANDS)
  character(len=256) :: srf_file = 's2a.srf'
  real(kind=8)       :: s2_response_mat(nw1nm, nb_S2)

contains

  subroutine semid_init(srf_fname)
    implicit none
    character(len=*), optional :: srf_fname
    if( present(srf_fname) ) then
       call set_srf_file(srf_fname)
    endif

    call read_s2_response()

  end subroutine semid_init

  subroutine set_srf_file(fname)
    implicit none
    character(len=*), intent(in) :: fname
    if( len(fname).gt.256 ) then
       write(*, '(a)') ' FATAL::set_srf_file:length of filename ***'//fname//'***'//&
            ' exceeds limit=256.'
       stop
    else if( len_trim(fname).eq.0 ) then
       write(*, '(a)') ' FATAL::set_srf_file:empty filename ***'//fname//'***'//&
            ' ?'
       stop
    else
       srf_file = fname
    endif
  end subroutine set_srf_file


  subroutine read_s2_response()
    implicit none
    integer :: rc
    rc = load_surface_data(srf_file, nb_S2, nw1nm, s2_response_mat)
    if( rc.ne.0 ) then
       write(*, '(a)' ) ' FATAL::read_s2_response: S2 response functions could not be '//&
            'read from file ***'//trim(srf_file)//'***'
       stop
    endif
  end subroutine read_s2_response


  !***********************************************************
  !     load_surface_datamodel_semidis_1geom
  !
  !> @brief read spectral surface response functions of Sensor from ASCII text file.
  !
  !> @details ordering is expected to be: Sensor bands column-wise, wavelenghts row-wise
  !
  !> @param[in]  fname     file name of ASCII text file providing spectral responses
  !> @param[in]  nCols     number of columns in file (expecting: number of Sensor bands)
  !> @param[in]  nRows     number of rows in file (expecting: wavelengths)
  !> @param[out] surf_mat  matrix of size nrows*ncols that will eventually hold the spectral responses
  !
  integer function load_surface_data(fname, nCols, nRows, surf_mat)
    implicit none
    ! arguments
    character(len=*), intent(in) :: fname
    integer, intent(in) :: nCols, nRows
    real(kind=8) :: surf_mat(nrows,ncols)
    ! local decls
    integer, parameter :: iounit = 18
    logical :: exists
    integer :: iostat
    integer :: nl,il

    !-- initialise
    load_surface_data = -1

    !-- check file
    inquire(FILE=fname, EXIST=exists)
    if( .not.exists ) then
       write(*, '(a)') ' FATAL::load_surface_data: ***'//trim(fname)//'*** was not found.'
    endif

    nl = cnt_lines(fname)
    if( nl.ne.nRows ) then
       write(*, '(2(a,i6))') ' ERROR::load_surface_data: detected nlines=',nl,&
            ' but expected nRows=',nRows
       return
    endif

    !-- open file
    open(iounit, file=fname)

    load_surface_data = 0
    row_loop:do il=1,nl
       read(iounit, *, iostat=iostat) surf_mat(il,:)
       if( iostat/=0 ) then
          write(*, '(a,i5)') ' FATAL::error occured when reading line=',il
          load_surface_data = -1
          exit row_loop
       endif
    enddo row_loop

    !-- close file
    close(iounit)
  contains
    integer function cnt_lines(fname)
      implicit none
      character(len=*), intent(in) :: fname
      integer, parameter :: lunit = 18

      cnt_lines = 0
      open(lunit, file=fname)
      do
         read (lunit,*, end=10)
         cnt_lines = cnt_lines + 1
      end do
10    close (lunit)
    end function cnt_lines

  end function load_surface_data

end module mo_sensimul_s2

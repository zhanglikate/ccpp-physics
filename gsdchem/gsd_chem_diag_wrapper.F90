!>\file gsd_chem_diag_wrapper.F90
!! This file is GSD Chemistry diagnostic wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 08/2020

 module gsd_chem_diag_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use gsd_chem_config
   use gocart_diag_mod

   implicit none

   private

   public :: gsd_chem_diag_wrapper_init, gsd_chem_diag_wrapper_run, gsd_chem_diag_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine gsd_chem_diag_wrapper_init()
      end subroutine gsd_chem_diag_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_gsd_chem_diag_wrapper_finalize Argument Table
!!
      subroutine gsd_chem_diag_wrapper_finalize()
      end subroutine gsd_chem_diag_wrapper_finalize

!> \defgroup gsd_chem_diag_group GSD Chem diag wrapper Module
!! This is the gsd chemistry
!>\defgroup gsd_chem_diag_wrapper GSD Chem diag wrapper Module  
!> \ingroup gsd_chem_diag_group
!! This is the GSD Chem diag wrapper Module
!! \section arg_table_gsd_chem_diag_wrapper_run Argument Table
!! \htmlinclude gsd_chem_diag_wrapper_run.html
!!
!>\section gsd_chem_diag_wrapper GSD Chemistry Scheme General Algorithm
!> @{
    subroutine gsd_chem_diag_wrapper_run(im, kte, kme, ktau,               &
                   pr3d, ntrac, ntso2, gq0,  tile_num, tmpmax, tmpmin,     &
                   spfhmax, spfhmin, t2m,   &
                   ntchmdiag, nseasalt,drydep, wetdpl, ssem, & 
                   errmsg,errflg)

    implicit none


    integer,        intent(in) :: im,kte,kme,ktau,tile_num
    integer,        intent(in) :: ntrac,ntso2,ntchmdiag, nseasalt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im), intent(inout) :: tmpmax, tmpmin, spfhmax, spfhmin, t2m
    real(kind_phys), dimension(im,kme), intent(in) :: pr3d
    real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: drydep, wetdpl
    real(kind_phys), dimension(im,nseasalt), intent(inout) :: ssem
    real(kind_phys), dimension(im,kte,ntrac), intent(in) :: gq0
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! -- for diagnostics
    real(kind_phys), dimension(ims:im, jms:jme, 6) :: trcm  ! inst tracer column mass density
    real(kind_phys), dimension(im,jme,kte,ntrac) :: gq0j
    real(kind_phys), dimension(im,jme,kme) :: pr3dj

    integer :: ide, ime, ite, kde


!>-- local variables
    integer :: i, j, jp, k, kp, nbegin
  

    errmsg = ''
    errflg = 0

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    nbegin = ntso2-1
    pr3dj(:,1,:  )=pr3d(:,:  )
    gq0j (:,1,:,:)=gq0 (:,:,:)
    ! -- calculate column mass density
    call gocart_diag_cmass(chem_opt, nbegin, g, pr3dj, gq0j, trcm)

    do i=1,im
     tmpmax(i)=drydep(i,p_seas_1) 
     !print*,'hli2 ssem',ssem(i,1),ssem(i,2),ssem(i,3),ssem(i,4),ssem(i,5)
     tmpmin(i)=(ssem(i,1)+ssem(i,2)+ssem(i,3)+ssem(i,4)+ssem(i,5))
     spfhmax(i)=wetdpl(i,p_seas_1) 
     spfhmin(i)=wetdpl(i,p_oc1) 
     t2m(i)=trcm(i,1,6) ! sea salt
    enddo
!   call gsd_chem_post() ! postprocessing for diagnostics

!
   end subroutine gsd_chem_diag_wrapper_run
!> @}

  end module gsd_chem_diag_wrapper

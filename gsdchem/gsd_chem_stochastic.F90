!>\file gsd_chem_seas_wrapper.F90
!! This file initializes and updates the stochastic perturbations array for emissions.
!! Unlike other gsd_chem-ccpp connection modules, this one is entirely standalone;
!! it uses no other gsd-chem modules.
!! Samuel.Trahan@noaa.gov 12/2020
module gsd_chem_stochastic
  use machine, only: kind_phys
  implicit none
  private

  public gsd_chem_stochastic_init, gsd_chem_stochastic_finalize, gsd_chem_stochastic_run

contains

  !> \brief Initialize the emissions multiplier to 1.0, indicating no
  !> emissions perturbations.
  !!
  subroutine gsd_chem_stochastic_init(im,emis_multiplier,errmsg,errflg)
    implicit none
    real(kind_phys), intent(out) :: emis_multiplier(im)
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    integer, intent(in) :: im

    errmsg=''
    errflg=0

    emis_multiplier = 1.0
  end subroutine gsd_chem_stochastic_init

  !> \brief Does absolutely nothing, and does it well.
  !!
  subroutine gsd_chem_stochastic_finalize
  end subroutine gsd_chem_stochastic_finalize

  !> \brief Update the emissions multiplier
  !! Updates the emis_multiplier with the output of
  !! stochastic_physics, if stochastic emissions perturbations are
  !! enabled.
  !!
  subroutine gsd_chem_stochastic_run(im, kme, emis_multiplier, ca1, ca_global_emis, &
                                     do_sppt_emis, sppt_wts, errmsg, errflg)
    implicit none

    logical,          intent(in) :: ca_global_emis, do_sppt_emis
    real(kind_phys),  intent(inout) :: emis_multiplier(im)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    integer,          intent(in) :: im,kme
    real(kind_phys), optional, intent(in) :: ca1(im)
    real(kind_phys), optional, intent(in) :: sppt_wts(:,:)

    integer :: i
    real(kind_phys) :: ca1_scaled

    if (do_sppt_emis) then
      do i = 1, im
        emis_multiplier(i) = max(0.5,min(1.5,sppt_wts(i,kme/2)))
      enddo
    elseif (ca_global_emis) then
      do i = 1, im
        ! ca1(i) is always precisely 0 or 2
        if(ca1(i)<1.0) then
          ca1_scaled=0.9
        else
          ca1_scaled=1.0/0.9
        endif
        emis_multiplier(i) = max(0.5,min(1.5,emis_multiplier(i)*0.95 + ca1_scaled*0.05))
      enddo
    endif
  end subroutine gsd_chem_stochastic_run

end module gsd_chem_stochastic

module mech_bridge
  ! Temporary bridge: wraps ddm_mechanics so we can later insert
  ! natural‐fracture effects without touching the drivers again.
  use ddm_mechanics, only: aperture_from_pressure_ddm => aperture_from_pressure, &
                           get_elem_stress_ddm => get_elem_stress
  use fracture_types
  implicit none
  private
  public :: mech_update, get_elem_stress
contains
  !---------------------------------------------------------------
  subroutine mech_update(muF, dt, ds)
    ! Main entry called by drivers:
    ! 1) Gather pressures of active HF elements (assumed to be the ones
    !    for which is_hf==true) preserving their original ordering so
    !    they match the existing DDM matrices.
    ! 2) Call the old DDM compliance to obtain apertures for those HF
    !    faces.
    ! 3) For active NF faces use a linear stiffness w = aperture0 + p/Kn
    ! 4) Push the updated apertures back to frac_list(:)%aperture.

    real*8, intent(in) :: muF, dt, ds   ! (placeholders – not yet used)

    if (.not.allocated(frac_list)) return

    integer :: nTot, nHF, i, h, nf
    nTot = size(frac_list)

    ! count HF active faces (they keep the original ordering)
    nHF = 0
    do i=1,nTot
       if (frac_list(i)%is_hf) nHF = nHF + 1
    end do

    if (nHF > 0) then
       real*8, allocatable :: p_hf(:), w_hf(:)
       allocate(p_hf(nHF), w_hf(nHF))

       h = 0
       do i=1,nTot
          if (frac_list(i)%is_hf) then
             h = h + 1
             p_hf(h) = frac_list(i)%pressure
          end if
       end do

       call aperture_from_pressure_ddm(p_hf, w_hf)

       h = 0
       do i=1,nTot
          if (frac_list(i)%is_hf) then
             h = h + 1
             frac_list(i)%aperture = w_hf(h)
          end if
       end do
    end if

    ! Natural fractures – linear elastic normal compliance
    do i=1,nTot
       if (.not.frac_list(i)%is_hf .and. frac_list(i)%is_active) then
          if (frac_list(i)%Kn > 0.0_dp) then
             frac_list(i)%aperture = frac_list(i)%aperture0 + &
                                     frac_list(i)%pressure / frac_list(i)%Kn
          else
             frac_list(i)%aperture = frac_list(i)%aperture0
          end if
       end if
    end do

  end subroutine mech_update

  !---------------------------------------------------------------
  ! Lightweight pass-through wrapper so other modules can query the
  ! element-centred stress tensor without importing ddm_mechanics
  ! directly (keeps coupling minimal).
  !---------------------------------------------------------------
  subroutine get_elem_stress(idx, sigma_cart)
    use fracture_types, only: dp
    implicit none
    integer, intent(in)   :: idx
    real(dp), intent(out) :: sigma_cart(3,3)
    call get_elem_stress_ddm(idx, sigma_cart)
  end subroutine get_elem_stress
end module mech_bridge 
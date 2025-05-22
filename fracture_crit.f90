module fracture_crit
  implicit none
contains
  !------------------------------------------------------------------
  ! Very simple tensile failure check (kept for backward compatibility)
  !------------------------------------------------------------------
  subroutine update_fracture_status(num_faces, sigma_normal, p, t0, status)
    integer, intent(in) :: num_faces
    real*8, intent(in) :: sigma_normal(num_faces)
    real*8, intent(in) :: p(num_faces)
    real*8, intent(in) :: t0
    integer, intent(inout) :: status(num_faces)
    integer :: i
    do i = 1, num_faces
      if (sigma_normal(i) + p(i) >= t0) then
        status(i) = 1
      end if
    end do
  end subroutine update_fracture_status

  !------------------------------------------------------------------
  ! Estimate mode-I Stress-Intensity Factor (K_I) at each face.
  ! A pragmatic approximation:  K_I = (p + σ_n) * sqrt(pi * w)
  ! where w = aperture.  This is NOT exact DDM LEFM but serves as a
  ! placeholder until the full singular integration is wired in.
  !------------------------------------------------------------------
  subroutine compute_KI(num_faces, pressure, sigma_normal, aperture, KI)
    integer, intent(in)  :: num_faces
    real*8 , intent(in)  :: pressure(num_faces)
    real*8 , intent(in)  :: sigma_normal(num_faces)
    real*8 , intent(in)  :: aperture(num_faces)
    real*8 , intent(out) :: KI(num_faces)
    real*8 , parameter   :: pi = 3.141592653589793d0
    integer :: i
    do i = 1, num_faces
       KI(i) = (pressure(i) + sigma_normal(i)) * sqrt(pi * max(aperture(i),1.0d-12))
    end do
  end subroutine compute_KI

  !------------------------------------------------------------------
  ! Update fracture status based on LEFM criterion K_I ≥ K_Ic
  !------------------------------------------------------------------
  subroutine update_fracture_lefm(num_faces, KI, KIc, status)
    integer, intent(in)  :: num_faces
    real*8 , intent(in)  :: KI(num_faces)
    real*8 , intent(in)  :: KIc
    integer, intent(inout) :: status(num_faces)
    integer :: i
    do i = 1, num_faces
       if (KI(i) >= KIc) status(i) = 1
    end do
  end subroutine update_fracture_lefm
end module fracture_crit

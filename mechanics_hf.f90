module mechanics_hf
  implicit none
contains
  subroutine update_mechanics(num_faces, pressure, krock, aperture, sigma_n)
    ! Simple linear normal-stiffness response for crack faces
    integer, intent(in) :: num_faces
    real*8, intent(in) :: pressure(num_faces)
    real*8, intent(in) :: krock
    real*8, intent(inout) :: aperture(num_faces)
    real*8, intent(out)  :: sigma_n(num_faces)
    integer :: i

    do i = 1, num_faces
       aperture(i) = pressure(i) / max(krock, 1.0d-12)
       sigma_n(i)  = krock * aperture(i)
    end do
  end subroutine update_mechanics
end module mechanics_hf

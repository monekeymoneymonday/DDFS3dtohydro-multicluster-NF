module hf_front
  !------------------------------------------------------------------
  ! Minimal front-tracking utility: grows the hydraulic fracture by
  ! cloning the current tip quadrilateral and translating it along the
  ! local x-axis (propagation direction) by a user-supplied distance ds.
  ! The new element is appended to frac_list and neighbour connectivity
  ! is rebuilt for involved faces.
  !------------------------------------------------------------------
  use fracture_types
  use fracture_utils, only: compute_neighbors
  use ki_ddm,        only: convert      ! local basis helper
  implicit none
  private
  public :: grow_hf_tip
contains

  subroutine grow_hf_tip(tip_idx, ds, new_idx)
    integer, intent(in)  :: tip_idx   ! index of current HF tip element in frac_list
    real(dp), intent(in) :: ds        ! propagation step length [m]
    integer, intent(out) :: new_idx   ! index of newly created element in frac_list

    real(dp) :: IEV(3,3), vec(3)
    real(dp) :: x(4), y(4), z(4)
    integer  :: i

    if (.not.allocated(frac_list)) stop 'grow_hf_tip: frac_list not allocated'
    if (tip_idx < 1 .or. tip_idx > size(frac_list)) stop 'grow_hf_tip: invalid tip idx'

    ! Build local coordinate system using first 3 vertices
    call convert( frac_list(tip_idx)%x(1), frac_list(tip_idx)%y(1), frac_list(tip_idx)%z(1), &
                 frac_list(tip_idx)%x(2), frac_list(tip_idx)%y(2), frac_list(tip_idx)%z(2), &
                 frac_list(tip_idx)%x(3), frac_list(tip_idx)%y(3), frac_list(tip_idx)%z(3), IEV )

    vec = ds * IEV(1,:)   ! local x-axis = propagation direction

    do i = 1, 4
       x(i) = frac_list(tip_idx)%x(i) + vec(1)
       y(i) = frac_list(tip_idx)%y(i) + vec(2)
       z(i) = frac_list(tip_idx)%z(i) + vec(3)
    end do

    ! Create new HF element via fracture_utils helper
    integer :: tmp_idx
    call create_hf_quad(x, y, z, frac_list(tip_idx)%aperture0, &
                        frac_list(tip_idx)%Kn, frac_list(tip_idx)%Ks, &
                        frac_list(tip_idx)%mu, frac_list(tip_idx)%cohesion, &
                        frac_list(tip_idx)%tensile_strength, .true., tmp_idx)

    new_idx = tmp_idx

    ! Recompute neighbours for the tip and new element only (cheap)
    call compute_neighbors()
  end subroutine grow_hf_tip

end module hf_front 
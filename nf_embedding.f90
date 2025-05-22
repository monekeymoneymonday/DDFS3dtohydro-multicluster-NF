module nf_embedding
  !---------------------------------------------------------------
  ! Utilities to find which fracture‐mesh element corresponds to
  ! each natural fracture from nf_list, based on coincident vertices.
  !---------------------------------------------------------------
  use fracture_types
  use nf_types, only: nf_list
  implicit none
  private
  public :: nf_elem_map, build_nf_elem_map

  integer, allocatable, public :: nf_elem_map(:) ! nf_elem_map(i)=element index in frac_list for NF i; 0 if not found

contains

  subroutine build_nf_elem_map(tol)
    real(dp), intent(in), optional :: tol
    real(dp) :: tol_use
    integer :: nNF, nEle, iNF, iEl, v

    if (.not.allocated(nf_list)) return
    if (.not.allocated(frac_list)) return

    tol_use = merge(tol, 1.0e-6_dp, present(tol))
    nNF  = size(nf_list)
    nEle = size(frac_list)

    if (.not.allocated(nf_elem_map)) allocate(nf_elem_map(nNF))
    nf_elem_map = 0

    do iNF = 1, nNF
       do iEl = 1, nEle
          if (.not.frac_list(iEl)%is_hf) cycle   ! only HF mesh for now
          if (vertices_coincident(iNF, iEl, tol_use)) then
             nf_elem_map(iNF) = iEl
             exit
          end if
       end do
    end do
  end subroutine build_nf_elem_map

  !-------------------------------------------------------------
  pure logical function vertices_coincident(idxNF, idxEl, tol)
    integer, intent(in) :: idxNF, idxEl
    real(dp), intent(in) :: tol
    logical :: used(4)
    integer :: i, j

    vertices_coincident = .false.
    used = .false.
    do i = 1, 4
       logical :: found
       found = .false.
       do j = 1, 4
          if (used(j)) cycle
          if (point_equal(nf_list(idxNF)%x(i), nf_list(idxNF)%y(i), nf_list(idxNF)%z(i), &
                          frac_list(idxEl)%x(j), frac_list(idxEl)%y(j), frac_list(idxEl)%z(j), tol)) then
             used(j) = .true.
             found = .true.
             exit
          end if
       end do
       if (.not.found) return
    end do
    vertices_coincident = .true.
  end function vertices_coincident

  pure logical function point_equal(x1,y1,z1,x2,y2,z2,tol)
    real(dp), intent(in):: x1,y1,z1,x2,y2,z2,tol
    point_equal = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 ) < tol
  end function point_equal

end module nf_embedding 
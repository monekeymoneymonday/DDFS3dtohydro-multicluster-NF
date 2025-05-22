module hf_nf_interaction
  !----------------------------------------------------------------
  ! Proximity detection between newly created HF elements and the
  ! catalogue of natural fractures. When an HF face is within a user
  ! tolerance of an NF plane, an event is pushed into a global queue
  ! for later evaluation (slip/open/cross logic).
  !----------------------------------------------------------------
  use fracture_types
  use nf_types,       only: nf_list
  implicit none
  private
  public :: InteractionEvent, event_queue, check_proximity

  type :: InteractionEvent
     integer :: hf_elem   = 0   ! index in frac_list
     integer :: nf_id     = 0   ! index in nf_list
     real(dp) :: distance = 0.0_dp
  end type InteractionEvent

  type(InteractionEvent), allocatable, public :: event_queue(:)

contains

  subroutine push_event(hf_idx, nf_idx, dist)
    integer, intent(in) :: hf_idx, nf_idx
    real(dp), intent(in) :: dist
    type(InteractionEvent) :: ev
    ev%hf_elem = hf_idx
    ev%nf_id   = nf_idx
    ev%distance = dist

    if (.not.allocated(event_queue)) then
       allocate(event_queue(1))
       event_queue(1) = ev
    else
       type(InteractionEvent), allocatable :: tmp(:)
       call move_alloc(event_queue, tmp)
       allocate(event_queue(size(tmp)+1))
       event_queue(1:size(tmp)) = tmp
       event_queue(size(tmp)+1) = ev
    end if
  end subroutine push_event

  !----------------------------------------------------------------
  subroutine check_proximity(new_hf_idx, d_tol)
    integer, intent(in) :: new_hf_idx
    real(dp), intent(in), optional :: d_tol
    real(dp) :: tol
    integer :: iNF

    if (.not.allocated(nf_list)) return
    if (.not.allocated(frac_list)) return

    tol = merge(d_tol, 1.0_dp, present(d_tol))

    ! centroid of HF element
    real(dp) :: pH(3)
    pH(1) = sum(frac_list(new_hf_idx)%x)/4.d0
    pH(2) = sum(frac_list(new_hf_idx)%y)/4.d0
    pH(3) = sum(frac_list(new_hf_idx)%z)/4.d0

    do iNF=1,size(nf_list)
       real(dp) :: p0(3), nvec(3), v1(3), v2(3), normN, dist
       p0 = [ sum(nf_list(iNF)%x)/4.d0, &
               sum(nf_list(iNF)%y)/4.d0, &
               sum(nf_list(iNF)%z)/4.d0 ]
       v1 = [ nf_list(iNF)%x(2)-nf_list(iNF)%x(1), &
               nf_list(iNF)%y(2)-nf_list(iNF)%y(1), &
               nf_list(iNF)%z(2)-nf_list(iNF)%z(1) ]
       v2 = [ nf_list(iNF)%x(3)-nf_list(iNF)%x(1), &
               nf_list(iNF)%y(3)-nf_list(iNF)%y(1), &
               nf_list(iNF)%z(3)-nf_list(iNF)%z(1) ]
       nvec(1) = v1(2)*v2(3) - v1(3)*v2(2)
       nvec(2) = v1(3)*v2(1) - v1(1)*v2(3)
       nvec(3) = v1(1)*v2(2) - v1(2)*v2(1)
       normN = sqrt(sum(nvec**2))
       if (normN==0.d0) cycle
       nvec = nvec / normN
       dist = abs( (pH(1)-p0(1))*nvec(1) + (pH(2)-p0(2))*nvec(2) + (pH(3)-p0(3))*nvec(3) )
       if (dist < tol) call push_event(new_hf_idx, iNF, dist)
    end do
  end subroutine check_proximity

end module hf_nf_interaction 
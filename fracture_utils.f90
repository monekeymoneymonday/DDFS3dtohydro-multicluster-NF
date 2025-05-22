module fracture_utils
  use fracture_types, only: dp, frac_elem, frac_list
  implicit none
contains
  !------------------------------------------------------------------
  ! Load centroid coordinates written by FSM3D_qua_part1 (output_part1.txt)
  ! The file has header lines followed by rows:
  !   i  Gx  Gy  Gz  Px  Py  Pz
  ! This routine extracts Gx, Gy, Gz for i = 1..num
  !------------------------------------------------------------------
  subroutine load_centers(filename, num, gx, gy, gz, ierr)
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: num
    real(dp),         intent(out) :: gx(num), gy(num), gz(num)
    integer,          intent(out) :: ierr

    character(len=256) :: line
    integer :: ios, idx, count
    real(dp)  :: x, y, z, dummy1, dummy2, dummy3

    ierr  = 0
    gx    = 0.0_dp
    gy    = 0.0_dp
    gz    = 0.0_dp
    count = 0

    open(unit=98,file=filename,status='old',action='read',iostat=ios)
    if (ios /= 0) then
       ierr = 1
       return
    end if

    do
       read(98,'(A)',iostat=ios) line
       if (ios /= 0) exit   ! EOF

       ! try to parse numeric line
       read(line,*,iostat=ios) idx, x, y, z, dummy1, dummy2, dummy3
       if (ios == 0) then
          if (idx >= 1 .and. idx <= num) then
             gx(idx) = x
             gy(idx) = y
             gz(idx) = z
             count   = count + 1
          end if
       end if
    end do
    close(98)

    if (count /= num) ierr = 2
  end subroutine load_centers

  !------------------------------------------------------------------
  ! Compute the average spacing between consecutive centroids
  !------------------------------------------------------------------
  function average_spacing(num, gx, gy, gz) result(ds)
    integer, intent(in) :: num
    real(dp), intent(in) :: gx(num), gy(num), gz(num)
    real(dp) :: ds, dist
    integer :: i

    if (num <= 1) then
       ds = 1.0_dp
       return
    end if

    ds = 0.0_dp
    do i = 1, num-1
       dist = sqrt( (gx(i+1)-gx(i))**2 + (gy(i+1)-gy(i))**2 + (gz(i+1)-gz(i))**2 )
       ds = ds + dist
    end do
    ds = ds / real(num-1, dp)
  end function average_spacing

  !================================================================
  ! PUBLIC API – element creation & connectivity
  !================================================================
  interface create_hf_quad
     module procedure create_quad_hf
  end interface
  interface create_nf_quad
     module procedure create_quad_nf
  end interface

  public :: create_hf_quad, create_nf_quad, append_element, compute_neighbors, update_leak_nf

  !----------------------------------------------------------------
  ! Append a fracture element to the global list
  !----------------------------------------------------------------
  subroutine append_element(elem, idx_out)
    type(frac_elem), intent(in)  :: elem
    integer,         intent(out) :: idx_out

    type(frac_elem), allocatable :: tmp(:)

    if (.not.allocated(frac_list)) then
       allocate(frac_list(1))
       frac_list(1) = elem
       idx_out = 1
       return
    end if

    call move_alloc(frac_list, tmp)           ! tmp now owns old data
    allocate(frac_list(size(tmp)+1))
    frac_list(1:size(tmp)) = tmp
    frac_list(size(tmp)+1) = elem
    idx_out = size(frac_list)
  end subroutine append_element

  !----------------------------------------------------------------
  ! Helper: compute unit normal from first three points of quadril.
  !----------------------------------------------------------------
  pure function compute_normal(x,y,z) result(n)
    real(dp), intent(in) :: x(4), y(4), z(4)
    real(dp) :: n(3), v1(3), v2(3), norm

    v1 = (/ x(2)-x(1), y(2)-y(1), z(2)-z(1) /)
    v2 = (/ x(3)-x(1), y(3)-y(1), z(3)-z(1) /)

    n(1) = v1(2)*v2(3) - v1(3)*v2(2)
    n(2) = v1(3)*v2(1) - v1(1)*v2(3)
    n(3) = v1(1)*v2(2) - v1(2)*v2(1)
    norm = sqrt(sum(n**2))
    if (norm > 0.0_dp) n = n / norm
  end function compute_normal

  !----------------------------------------------------------------
  ! Generic quad‐element creator used by HF & NF wrappers
  !----------------------------------------------------------------
  subroutine create_quad_generic(x,y,z, is_hf, aperture0, Kn, Ks, mu, cohesion, &
                                 tensile_strength, can_grow, idx_out)
    real(dp), intent(in) :: x(4), y(4), z(4)
    logical , intent(in) :: is_hf, can_grow
    real(dp), intent(in), optional :: aperture0, Kn, Ks, mu, cohesion, tensile_strength
    integer , intent(out) :: idx_out

    type(frac_elem) :: e

    e%x = x; e%y = y; e%z = z
    e%normal = compute_normal(x,y,z)

    if (present(aperture0))        e%aperture0 = aperture0
    if (present(Kn))               e%Kn        = Kn
    if (present(Ks))               e%Ks        = Ks
    if (present(mu))               e%mu        = mu
    if (present(cohesion))         e%cohesion  = cohesion
    if (present(tensile_strength)) e%tensile_strength = tensile_strength

    e%aperture = e%aperture0
    e%is_hf    = is_hf
    e%can_grow = can_grow
    e%is_active= is_hf             ! HF faces start active by default; NF inactive

    idx_out = -1
    call append_element(e, idx_out)
  end subroutine create_quad_generic

  ! HF wrapper
  subroutine create_quad_hf(x,y,z, aperture0, Kn, Ks, mu, cohesion, tensile_strength, &
                            can_grow, idx)
    real(dp), intent(in) :: x(4), y(4), z(4)
    real(dp), intent(in), optional :: aperture0, Kn, Ks, mu, cohesion, tensile_strength
    logical , intent(in), optional :: can_grow
    integer , intent(out) :: idx

    logical :: grow
    if (present(can_grow)) then
       grow = can_grow
    else
       grow = .false.
    end if
    call create_quad_generic(x,y,z,.true., aperture0, Kn, Ks, mu, cohesion, &
                             tensile_strength, grow, idx)
  end subroutine create_quad_hf

  ! NF wrapper
  subroutine create_quad_nf(x,y,z, aperture0, Kn, Ks, mu, cohesion, tensile_strength, idx)
    real(dp), intent(in) :: x(4), y(4), z(4)
    real(dp), intent(in), optional :: aperture0, Kn, Ks, mu, cohesion, tensile_strength
    integer , intent(out) :: idx

    call create_quad_generic(x,y,z,.false., aperture0, Kn, Ks, mu, cohesion, &
                             tensile_strength, .false., idx)
  end subroutine create_quad_nf

  !----------------------------------------------------------------
  ! Compute neighbour connectivity for all elements in frac_list.
  ! Neighbours share an edge (two common vertices within tolerance).
  !----------------------------------------------------------------
  subroutine compute_neighbors(tol)
    real(dp), intent(in), optional :: tol
    real(dp) :: tol_use
    integer :: i,j

    if (present(tol)) then
       tol_use = tol
    else
       tol_use = 1.0e-6_dp
    end if

    if (.not.allocated(frac_list)) return

    ! First clear any existing neighbour info
    do i = 1, size(frac_list)
       if (allocated(frac_list(i)%neighbors)) deallocate(frac_list(i)%neighbors)
    end do

    do i = 1, size(frac_list)-1
       do j = i+1, size(frac_list)
          if (share_edge(i,j,tol_use)) then
             call add_neighbor_to(i, j)
             call add_neighbor_to(j, i)
          end if
       end do
    end do
  end subroutine compute_neighbors

  !----------------------------------------------------------------
  ! Internal: test if two elements share an edge
  !----------------------------------------------------------------
  logical function share_edge(i,j,tol)
    integer, intent(in) :: i,j
    real(dp), intent(in) :: tol
    logical :: edge_match
    integer :: ei,ej, ei_nxt, ej_nxt
    real(dp) :: xi(4), yi(4), zi(4), xj(4), yj(4), zj(4)
    real(dp) :: p1(3), p2(3), q1(3), q2(3)

    xi = frac_list(i)%x; yi = frac_list(i)%y; zi = frac_list(i)%z
    xj = frac_list(j)%x; yj = frac_list(j)%y; zj = frac_list(j)%z

    share_edge = .false.

    do ei = 1,4
       ei_nxt = mod(ei,4)+1      ! wrap 4→1
       p1 = (/ xi(ei), yi(ei), zi(ei) /)
       p2 = (/ xi(ei_nxt), yi(ei_nxt), zi(ei_nxt) /)

       do ej = 1,4
          ej_nxt = mod(ej,4)+1
          q1 = (/ xj(ej), yj(ej), zj(ej) /)
          q2 = (/ xj(ej_nxt), yj(ej_nxt), zj(ej_nxt) /)

          edge_match = (points_equal(p1,q1,tol) .and. points_equal(p2,q2,tol)) .or. &
                       (points_equal(p1,q2,tol) .and. points_equal(p2,q1,tol))
          if (edge_match) then
             share_edge = .true.
             return
          end if
       end do
    end do
  end function share_edge

  !----------------------------------------------------------------
  pure logical function points_equal(p,q,tol)
    real(dp), intent(in) :: p(3), q(3), tol
    points_equal = sqrt( (p(1)-q(1))**2 + (p(2)-q(2))**2 + (p(3)-q(3))**2 ) < tol
  end function points_equal

  !----------------------------------------------------------------
  subroutine add_neighbor_to(idx, neigh)
    integer, intent(in) :: idx, neigh
    integer, allocatable :: tmp(:)

    if (allocated(frac_list(idx)%neighbors)) then
       call move_alloc(frac_list(idx)%neighbors, tmp)
       allocate(frac_list(idx)%neighbors(size(tmp)+1))
       frac_list(idx)%neighbors(1:size(tmp)) = tmp
       frac_list(idx)%neighbors(size(tmp)+1) = neigh
    else
       allocate(frac_list(idx)%neighbors(1))
       frac_list(idx)%neighbors(1) = neigh
    end if
  end subroutine add_neighbor_to

  !--------------------------------------------------------------
  ! Update Carter leak-off rate for active NF elements.
  !   C_L – Carter coefficient [m/s^0.5]
  !   t   – current simulation time [s]
  !   dt  – time step [s]  (used to avoid division by zero)
  !--------------------------------------------------------------
  subroutine update_leak_nf(C_L, t, dt)
    real(dp), intent(in) :: C_L, t, dt
    integer :: i

    if (.not.allocated(frac_list)) return

    do i = 1, size(frac_list)
       if (frac_list(i)%is_active .and. .not.frac_list(i)%is_hf) then
          if (frac_list(i)%t_open < 0.0_dp) then
             frac_list(i)%t_open = t - dt   ! mark opening at previous step
          end if
          if (frac_list(i)%t_open > 0.0_dp) then
             frac_list(i)%leak = C_L / sqrt(max(t - frac_list(i)%t_open, dt))
          else
             frac_list(i)%leak = 0.0_dp
          end if
       end if
    end do
  end subroutine update_leak_nf

end module fracture_utils 
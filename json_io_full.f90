module json_io_full
  ! Requires the JSON-Fortran library (https://github.com/jacobwilliams/json-fortran)
  use json_module
  use fracture_utils, only: create_nf_quad, compute_neighbors
  use fracture_types, only: dp
  implicit none
  private
  public :: read_config, write_state, load_natural_fractures
contains
  !------------------------------------------------------------------
  subroutine read_config(file, NELE, mu, dt, nstep, inlet_type, inlet_val, t0, KIc, C_L)
    character(len=*), intent(in) :: file
    integer,          intent(out):: NELE, nstep, inlet_type
    real*8,           intent(out):: mu, dt, inlet_val, t0, KIc, C_L
    type(json_file) :: jf
    call json_file_init(jf, file, new=.false.)
    call json_get(jf, '/NELE',  NELE)
    call json_get(jf, '/mu',    mu)
    call json_get(jf, '/dt',    dt)
    call json_get(jf, '/nstep', nstep)
    call json_get(jf, '/inlet_type', inlet_type)
    call json_get(jf, '/inlet_val',  inlet_val)
    call json_get(jf, '/C_L',        C_L)
    call json_get(jf, '/t0',    t0)
    call json_get(jf, '/KIc',   KIc)
    call json_file_destroy(jf)

    !------------------------------------------------------------
    ! Additionally load natural fractures (if present) and build
    ! connectivity. This side-effect populates the global frac_list.
    !------------------------------------------------------------
    call load_natural_fractures(file)
  end subroutine read_config

  !------------------------------------------------------------------
  subroutine write_state(file, step, pressure, aperture)
    character(len=*), intent(in) :: file
    integer,          intent(in) :: step
    real*8,           intent(in) :: pressure(:), aperture(:)
    type(json_file) :: jf
    call json_file_init(jf, file, overwrite=.true.)
    call json_add(jf,'step', step)
    call json_add(jf,'pressure', pressure)
    call json_add(jf,'aperture', aperture)
    call json_write(jf)
    call json_file_destroy(jf)
  end subroutine write_state

  !------------------------------------------------------------------
  ! Read natural fractures from JSON file and append them to frac_list
  ! Expected JSON structure:
  !   "natural_fractures" : [
  !        { "nodes" : [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3],[x4,y4,z4]],
  !          "aperture":1e-4, "Kn":1e11, "Ks":1e11,
  !          "mu":0.6, "cohesion":5e6, "tensile_strength":2e6 }
  !        , ...
  !   ]
  ! Additional convenience: if "centroid", "strike", "dip", "half_length", "half_height"
  ! are provided instead of nodes, a rectangular quad is generated.
  !------------------------------------------------------------------
  subroutine load_natural_fractures(file)
    character(len=*), intent(in) :: file

    type(json_file) :: jf
    integer :: nNF, i, ierr

    call json_file_init(jf, file, new=.false.)

    ! How many entries?
    call json_get(jf, '/natural_fractures', n_elements=nNF, found=ierr)
    if (ierr /= 0 .or. nNF <= 0) then
       call json_file_destroy(jf)
       return  ! nothing to do
    end if

    do i = 1, nNF
       call parse_single_nf(jf, i)
    end do

    call json_file_destroy(jf)

    ! update connectivity once all elements are in place
    call compute_neighbors()
  end subroutine load_natural_fractures

  !----------------------------------------------------------------
  subroutine parse_single_nf(jf, idx)
    type(json_file), intent(in) :: jf
    integer,        intent(in) :: idx

    character(len=:), allocatable :: base
    base = '/natural_fractures('//trim(adjustl(int_to_str(idx)))//')'

    ! Check if "nodes" array exists using json_get found flag
    logical :: has_nodes
    integer :: ierr_test, dummyN
    call json_get(jf, base//'/nodes', n_elements=dummyN, found=ierr_test)
    has_nodes = (ierr_test == 0)

    real(dp) :: x(4), y(4), z(4)
    real(dp) :: aperture0, Kn, Ks, mu, cohesion, T0
    aperture0 = 1.0e-4_dp; Kn = 1.0e12_dp; Ks = 1.0e12_dp
    mu = 0.6_dp; cohesion = 0.0_dp; T0 = 0.0_dp

    if (has_nodes) then
       real(dp), allocatable :: nodes(:,:)
       call json_get(jf, base//'/nodes', nodes)
       if (size(nodes,1) /= 4 .or. size(nodes,2) /= 3) then
          stop 'parse_single_nf: nodes array must be 4x3'
       end if
       x = nodes(:,1); y = nodes(:,2); z = nodes(:,3)
    else
       ! Alternative: centroid + orientation + sizes -> build rectangle
       real(dp) :: cen(3), strike, dip, hl, hh
       call json_get(jf, base//'/centroid', cen)
       call json_get(jf, base//'/strike',  strike)
       call json_get(jf, base//'/dip',     dip)
       call json_get(jf, base//'/half_length', hl)
       call json_get(jf, base//'/half_height', hh)

       call rectangle_from_strike_dip(cen, strike, dip, hl, hh, x, y, z)
    end if

    call json_get_optional(jf, base//'/aperture', aperture0)
    call json_get_optional(jf, base//'/Kn',       Kn)
    call json_get_optional(jf, base//'/Ks',       Ks)
    call json_get_optional(jf, base//'/mu',       mu)
    call json_get_optional(jf, base//'/cohesion', cohesion)
    call json_get_optional(jf, base//'/tensile_strength', T0)

    integer :: new_idx
    call create_nf_quad(x,y,z, aperture0, Kn, Ks, mu, cohesion, T0, new_idx)
  end subroutine parse_single_nf

  !----------------------------------------------------------------
  ! Build rectangle corner points from centroid + orientation
  !----------------------------------------------------------------
  subroutine rectangle_from_strike_dip(cen, strike, dip, hl, hh, x, y, z)
    real(dp), intent(in)  :: cen(3), strike, dip, hl, hh
    real(dp), intent(out) :: x(4), y(4), z(4)

    real(dp), parameter :: deg2rad = 3.14159265358979323846_dp/180.0_dp
    real(dp) :: st, cp, sp, cd, sd, e1(3), e2(3), vL(3), vH(3), p0(3), p(4,3)

    st = strike*deg2rad
    cp = cos(st); sp = sin(st)
    cd = cos(dip*deg2rad); sd = sin(dip*deg2rad)

    ! Local axes: e1 along strike (horizontal), e2 dip direction (down)
    e1 = (/ sp, cp, 0.0_dp /)
    e2 = (/ -cp*sd, sp*sd, cd /)   ! points downward normal to e1 x e3

    vL = hl*e1
    vH = hh*e2

    p0 = cen - vL - vH
    p(1,:) = p0
    p(2,:) = cen + vL - vH
    p(3,:) = cen + vL + vH
    p(4,:) = cen - vL + vH

    x = p(:,1); y = p(:,2); z = p(:,3)
  end subroutine rectangle_from_strike_dip

  !----------------------------------------------------------------
  ! JSON optional helper (overloaded for real*8)
  !----------------------------------------------------------------
  subroutine json_get_optional(jf, path, val)
    type(json_file), intent(in) :: jf
    character(len=*), intent(in):: path
    real(dp), intent(inout)     :: val
    real(dp) :: tmp
    integer :: ierr
    call json_get(jf, path, tmp, found=ierr)
    if (ierr == 0) val = tmp   ! only update if present
  end subroutine json_get_optional

  ! integer to string helper
  pure function int_to_str(i) result(res)
    integer, intent(in) :: i
    character(len=20)   :: res
    write(res,'(I0)') i
  end function int_to_str
end module json_io_full 
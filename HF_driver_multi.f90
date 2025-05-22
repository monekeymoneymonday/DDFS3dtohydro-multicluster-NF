!*********************************************************************
! HF_driver_multi.f90  – multi-cluster hydraulic-fracture driver
!   * Uses the flow_hf_multi wrapper to solve Poiseuille flow separately
!     on each fracture cluster (contiguous element range).
!   * Reads cluster layout and inlet specification from a lightweight
!     text file (see sample at end of this source).
!   * Keeps the same mechanical & fracture-propagation workflow already
!     implemented in the single-cluster driver.
!*********************************************************************
program HF_driver_multi
  use flow_hf_multi        ! NEW – multi-cluster pressure solver
  use flow_network,       only: solve_hf_pressure
  use mech_bridge          ! aperture from DDM compliance bridge
  use fracture_crit        ! failure criterion
  use fracture_utils       ! centroid loader & spacing calc (renamed)
  use ki_ddm               ! DDM-based KI computation
  use cluster_data         ! cluster parser (nClust, ranges, inlet arrays)
  use json_io              ! fallback lightweight dump
  use json_io_full, only : write_state
  use fracture_types       ! global fracture container
  implicit none

  !--------------- config / cluster data ---------------------------
  integer :: nClust                       ! number of clusters
  integer, allocatable :: iStart(:), iEnd(:)
  integer, allocatable :: inletType(:)
  real*8 , allocatable :: inletVal(:)
  integer :: NELE                         ! total number of faces

  !--------------- global control parameters -----------------------
  integer :: nstep
  real*8  :: mu_fluid, dt, t0, KIc, C_L

  !--------------- derived geometric data -------------------------
  real*8 :: ds
  real*8, allocatable :: gx(:), gy(:), gz(:)
  integer :: ierr, k, idx_init

  !--------------- primary fields ---------------------------------
  real*8 , allocatable :: aperture(:), pressure(:), sigma_n(:)
  integer, allocatable :: status(:), prev_status(:)
  real*8 , allocatable :: KI(:), t_open(:), leak(:)

  !--------------- misc -------------------------------------------
  integer :: step, pic, ic          ! loop over clusters
  real*8, parameter :: w0 = 1.0d-3   ! base aperture for broken faces [m]
  integer, parameter :: maxPic = 10
  real*8 , parameter :: tolPic = 1.0d-4

  !--------------- 0. Read configuration --------------------------
  call read_cluster_config('input_hf_multi.txt', nClust, iStart, iEnd,  &
                           inletType, inletVal, NELE)

  ! For now the remaining (scalar) controls are read from an auxiliary
  ! file named "hf_params.txt" with a fixed order, to keep the stub
  ! simple. The format is:
  !   mu_fluid  dt  nstep  t0  KIc  C_L
  open(unit=11,file='hf_params.txt',status='old',action='read')
  read(11,*) mu_fluid, dt, nstep, t0, KIc, C_L
  close(11)

  !--------------- 1. Average spacing from part-1 output ----------
  allocate(gx(NELE), gy(NELE), gz(NELE))
  call load_centers('../part1/output_part1.txt', NELE, gx, gy, gz, ierr)
  if (ierr == 0) then
     ds = average_spacing(NELE, gx, gy, gz)
  else
     ds = 1.0d0
  end if
  deallocate(gx, gy, gz)

  !--------------- 2. Allocate state arrays -----------------------
  allocate(aperture(NELE), pressure(NELE), sigma_n(NELE))
  allocate(status(NELE), prev_status(NELE))
  allocate(KI(NELE))
  allocate(t_open(NELE)); t_open = -1.0d0
  allocate(leak(NELE));   leak   = 0.0d0

  aperture = 1.0d-3
  pressure = 0.0d0
  sigma_n  = 0.0d0
  status   = 0
  prev_status = 0

  ! Initialise frac_list with HF faces (contiguous ordering)
  if (.not.allocated(frac_list)) allocate(frac_list(NELE))
  do idx_init = 1, NELE
     frac_list(idx_init)%is_hf = .true.
     frac_list(idx_init)%is_active = .true.
     frac_list(idx_init)%aperture0 = aperture(idx_init)
     frac_list(idx_init)%aperture  = aperture(idx_init)
     frac_list(idx_init)%pressure  = pressure(idx_init)
  end do

  !--------------- 3. Time-stepping loop --------------------------
  do step = 1, nstep

     ! (a) Mechanics update using last-known pressure
     do k = 1, NELE
        frac_list(k)%pressure = pressure(k)
     end do
     call mech_update(mu_fluid, dt, ds)
     do k = 1, NELE
        aperture(k) = frac_list(k)%aperture
     end do
     sigma_n = -pressure
     do k = 1, NELE
        if (status(k) == 1) aperture(k) = max(aperture(k), w0)
     end do

     ! (b) Carter leak-off
     call update_leakoff(step, dt, C_L, status, t_open, leak)

     ! (c) Picard coupling between mechanics & fluid ---------------
     call picard_multi(nClust, iStart, iEnd, inletType, inletVal,       &
                       aperture, mu_fluid, dt, ds, leak, pressure)

     ! (d) Fracture criterion & status update -----------------------
     KI = 0.0d0                             ! reset before fresh evaluation
     do ic = 1, nClust
        call get_KI_tip_ddm(iEnd(ic), KI(iEnd(ic)))
     end do
     call update_fracture_lefm(NELE, KI, KIc, status)

     ! mark newly broken faces traction-free
     do k = 1, NELE
        if (status(k)==1 .and. prev_status(k)==0) call set_face_broken(k)
     end do
     prev_status = status

     !------------------------------------------------------------
     ! Rolling JSON output (state.json) overwritten each step.
     ! Falls back to dump_state if JSON-Fortran unavailable.
     !------------------------------------------------------------
     logical, save :: use_full = .true.
     if (use_full) then
        call write_state('state.json', step, pressure, aperture)
     else
        call dump_state(step, pressure, aperture)
     end if

     !--------------------------------------------------------------
     ! Check every cluster tip; insert a new element right after the
     ! tip if it failed during this step. Order is preserved so the
     ! ranges remain contiguous.
     !--------------------------------------------------------------
     integer :: tip, oldN, j
     do ic = 1, nClust
        tip = iEnd(ic)
        if ( status(tip) == 1 ) then
           oldN = NELE
           !--- insert slot at position tip -----------------------
           call insert_real(aperture, tip, oldN)
           call insert_real(pressure, tip, oldN)
           call insert_real(sigma_n,  tip, oldN)
           call insert_int (status,   tip, oldN)
           call insert_int (prev_status, tip, oldN)
           call insert_real(KI,       tip, oldN)
           call insert_real(t_open,   tip, oldN)
           call insert_real(leak,     tip, oldN)

           NELE = oldN + 1

           ! initialise new face with tip values
           aperture(tip+1) = aperture(tip)
           pressure(tip+1) = pressure(tip)
           sigma_n(tip+1)  = sigma_n(tip)
           status(tip+1)   = 0
           prev_status(tip+1) = 0
           KI(tip+1)       = 0.0d0
           t_open(tip+1)   = -1.0d0
           leak(tip+1)     = 0.0d0

           ! Extend frac_list with the new HF element
           if (allocated(frac_list)) then
              type(frac_elem), allocatable :: tmpFL(:)
              call move_alloc(frac_list, tmpFL)
              allocate(frac_list(oldN+1))
              frac_list(1:oldN) = tmpFL
           else
              allocate(frac_list(NELE))
           end if
           frac_list(tip+1)%is_hf = .true.
           frac_list(tip+1)%is_active = .true.
           frac_list(tip+1)%aperture0 = aperture(tip+1)
           frac_list(tip+1)%aperture  = aperture(tip+1)
           frac_list(tip+1)%pressure  = pressure(tip+1)

           ! update cluster indices
           iEnd(ic) = iEnd(ic) + 1
           do j = ic+1, nClust
              iStart(j) = iStart(j) + 1
              iEnd(j)   = iEnd(j)   + 1
           end do

           ! create geometry (placeholder)
           call extend_mesh_insert(tip, ds)
        end if
     end do
  end do

  write(*,*) 'HF_driver_multi finished. Pressure field saved to pressure_out.txt'

contains
  !----------------------------------------------------------------
  subroutine update_leakoff(step, dt, C_L, status, t_open, leak)
    integer, intent(in) :: step
    real*8 , intent(in) :: dt, C_L
    integer, intent(in) :: status(:)
    real*8 , intent(inout) :: t_open(:), leak(:)
    real*8 :: current_time
    integer :: i
    current_time = step*dt
    do i = 1, size(status)
       if (status(i)==1 .and. t_open(i) < 0.0d0) t_open(i)=current_time-dt
       if (t_open(i) > 0.0d0) then
          leak(i) = C_L / sqrt(max(current_time - t_open(i), dt))
       else
          leak(i) = 0.0d0
       end if
    end do
  end subroutine update_leakoff

  !----------------------------------------------------------------
  subroutine picard_multi(nClust, iStart, iEnd, inletType, inletVal,    &
                          aperture, mu, dt, ds, leak, pressure)
    ! Iterate between mechanics and multi-cluster fluid flow until
    ! pressures converge for the whole domain.  For clusters marked
    ! inletType=1 (rate-controlled) an inner Newton loop adjusts the
    ! inlet pressure such that the prescribed volumetric rate is
    ! satisfied.

    integer, intent(in)    :: nClust, iStart(nClust), iEnd(nClust)
    integer, intent(in)    :: inletType(nClust)
    real*8 , intent(in)    :: inletVal(nClust)
    real*8 , intent(inout) :: aperture(:)
    real*8 , intent(in)    :: mu, dt, ds
    real*8 , intent(in)    :: leak(:)
    real*8 , intent(inout) :: pressure(:)

    real*8 , allocatable :: p_iter(:), p_new(:), p_tmp(:)
    integer :: pic, ic, i1, i2, nLoc, k

    if (.not.allocated(p_iter)) allocate(p_iter(size(pressure)), p_new(size(pressure)))
    p_iter = pressure

    do pic = 1, maxPic
       !-----------------------------------------------------------
       ! (1) Mechanics with current pressure guess
       !-----------------------------------------------------------
       do k = 1, size(pressure)
          frac_list(k)%pressure = p_iter(k)
       end do
       call mech_update(mu, dt, ds)
       do k = 1, size(pressure)
          aperture(k) = frac_list(k)%aperture
       end do

       !-----------------------------------------------------------
       ! (2) Fluid solve cluster-by-cluster (with rate control)
       !-----------------------------------------------------------
       p_new = p_iter   ! start from previous for safety

       do ic = 1, nClust
          i1   = iStart(ic)
          i2   = iEnd(ic)
          nLoc = i2 - i1 + 1
          if (nLoc <= 0) cycle

          if (inletType(ic) == 0) then
             ! ---- Dirichlet inlet pressure ----
             call solve_hf_pressure(nLoc, aperture(i1:i2), mu, dt, ds,      &
                                   p_iter(i1:i2), 0, inletVal(ic),          &
                                   leak(i1:i2), p_new(i1:i2))
          else
             ! ---- Rate-controlled inlet (Newton iteration on Pin) ----
             real*8 :: Pin, Qcalc, k12, Qtol
             integer, parameter :: maxIt = 10
             integer :: it
             logical :: newton_ok

             if (.not.allocated(p_tmp)) allocate(p_tmp(nLoc))

             !--- Guard against clusters shorter than 2 faces ---------
             if (nLoc < 2) then
                write(*,*) 'Warning: Rate-controlled cluster of length ', nLoc, &
                           ' - reverting to Dirichlet treatment.'
                call solve_hf_pressure(nLoc, aperture(i1:i2), mu, dt, ds,   &
                                      p_iter(i1:i2), 0, inletVal(ic),       &
                                      leak(i1:i2), p_new(i1:i2))
                cycle
             end if

             ! crude initial guess – use previous first-node pressure if available
             Pin = max(p_iter(i1), 1.0d3)
             Qtol = 1.d-6*inletVal(ic) + 1.d-12
             newton_ok = .false.

             do it = 1, maxIt
                call solve_hf_pressure(nLoc, aperture(i1:i2), mu, dt, ds,   &
                                      p_iter(i1:i2), 0, Pin,                &
                                      leak(i1:i2), p_tmp)

                ! permeability between first two faces
                k12   = ((aperture(i1) + aperture(i1+1)) / 2.d0)**3 / (12.d0*mu)
                Qcalc = k12 * (Pin - p_tmp(2)) / ds

                if (abs(Qcalc - inletVal(ic)) < Qtol) then
                   newton_ok = .true.
                   exit
                end if

                ! Newton update on inlet pressure
                Pin = Pin + (inletVal(ic) - Qcalc)*ds / k12
             end do

             if (.not.newton_ok) then
                write(*,*) 'Warning: Newton did not converge for cluster ', ic, &
                           ' after ', maxIt, ' iterations. Proceeding with last iterate.'
             end if

             p_new(i1:i2) = p_tmp
          end if
       end do   ! cluster loop

       !-----------------------------------------------------------
       ! (3) Convergence check over the whole domain
       !-----------------------------------------------------------
       if (maxval(abs(p_new - p_iter)) < tolPic*maxval(abs(p_iter)+1.d-12)) exit
       p_iter = p_new
    end do   ! Picard outer loop

    pressure = p_new
  end subroutine picard_multi

  !----------------------------------------------------------------
  ! Helper generics to append one slot to an allocatable array
  !----------------------------------------------------------------
  subroutine extend_real(arr, oldN)
    real*8, allocatable, intent(inout) :: arr(:)
    integer,            intent(in)    :: oldN
    real*8, allocatable :: tmp(:)
    call move_alloc(arr, tmp)
    allocate(arr(oldN+1))
    arr(1:oldN) = tmp
    arr(oldN+1) = 0.0d0
  end subroutine extend_real

  subroutine extend_int(arr, oldN)
    integer, allocatable, intent(inout) :: arr(:)
    integer,            intent(in)     :: oldN
    integer, allocatable :: tmp(:)
    call move_alloc(arr, tmp)
    allocate(arr(oldN+1))
    arr(1:oldN) = tmp
    arr(oldN+1) = 0
  end subroutine extend_int

  subroutine insert_real(arr, pos, oldN)
    real*8, allocatable, intent(inout) :: arr(:)
    integer,            intent(in)    :: pos, oldN
    real*8, allocatable :: tmp(:)
    call move_alloc(arr, tmp)
    allocate(arr(oldN+1))
    arr(1:pos) = tmp(1:pos)
    arr(pos+1:oldN+1) = tmp(pos:oldN)
  end subroutine insert_real

  subroutine insert_int(arr, pos, oldN)
    integer, allocatable, intent(inout) :: arr(:)
    integer,            intent(in)     :: pos, oldN
    integer, allocatable :: tmp(:)
    call move_alloc(arr, tmp)
    allocate(arr(oldN+1))
    arr(1:pos) = tmp(1:pos)
    arr(pos+1:oldN+1) = tmp(pos:oldN)
  end subroutine insert_int
end program HF_driver_multi 
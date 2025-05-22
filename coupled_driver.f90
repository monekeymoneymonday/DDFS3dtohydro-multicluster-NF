program coupled_driver
  !------------------------------------------------------------------
  ! Generic coupled HF + NF driver.
  !  – Supports several HF clusters (contiguous element ranges) that
  !    may interact with natural fractures.
  !  – Alternates mechanics (ddm compliance) and fluid flow until
  !    pressures converge (Picard iteration).
  !  – Carter leak-off applied both to HF and active NF faces.
  !------------------------------------------------------------------
  use fracture_types
  use fracture_utils,  only: load_centers, average_spacing, update_leak_nf
  use cluster_data      ! provides cluster layout & inlet settings
  use mech_bridge       ! aperture_from_pressure + stress transfer
  use flow_network      ! global network flow solver (handles NF leak)
  use fracture_crit     ! LEFM criterion (placeholder)
  use ki_ddm            ! KI computations
  use json_io_full, only: write_state
  implicit none

  ! ---------------- configuration ----------------------------------
  integer :: nClust
  integer, allocatable :: iStart(:), iEnd(:)       ! cluster ranges (1-based)
  integer, allocatable :: inletType(:)             ! 0=P fixed, 1=rate Q
  real*8 , allocatable :: inletVal(:)
  integer :: NELE

  ! Scalar params read from file
  integer :: nstep
  real*8  :: mu_fluid, dt, t0, KIc, C_L

  ! ---------------- geometry / spacing -----------------------------
  real*8 :: ds
  real*8, allocatable :: gx(:), gy(:), gz(:)
  integer :: ierr

  ! ---------------- state arrays -----------------------------------
  real*8 , allocatable :: aperture(:), pressure(:)
  real*8 , allocatable :: sigma_n(:)
  integer, allocatable :: status(:), prev_status(:)
  real*8 , allocatable :: KI(:)

  ! ---------------- misc -------------------------------------------
  integer :: step, pic, ic, k
  integer, parameter :: maxPic = 10
  real*8 , parameter :: tolPic = 1.0d-4
  real*8 , parameter :: w0     = 1.0d-3   ! base open width [m]

  ! ================================================================
  ! 1. Input --------------------------------------------------------
  ! ================================================================
  call read_cluster_config('input_hf_multi.txt', nClust, iStart, iEnd, &
                           inletType, inletVal, NELE)

  open(11,file='hf_params.txt',status='old',action='read')
  read(11,*) mu_fluid, dt, nstep, t0, KIc, C_L
  close(11)

  ! spacing estimate
  allocate(gx(NELE), gy(NELE), gz(NELE))
  call load_centers('../part1/output_part1.txt', NELE, gx, gy, gz, ierr)
  if (ierr==0) then
     ds = average_spacing(NELE, gx, gy, gz)
  else
     ds = 1.0d0
  end if
  deallocate(gx,gy,gz)

  ! ================================================================
  ! 2. Allocate arrays & initialise frac_list ----------------------
  ! ================================================================
  allocate(aperture(NELE), pressure(NELE), sigma_n(NELE))
  allocate(status(NELE), prev_status(NELE))
  allocate(KI(NELE))

  aperture = 1.0d-3; pressure = 0.0d0; sigma_n = 0.0d0
  status   = 0; prev_status = 0; KI = 0.0d0

  if (.not.allocated(frac_list)) allocate(frac_list(NELE))
  do k=1,NELE
     frac_list(k)%is_hf     = .true.
     frac_list(k)%is_active = .true.
     frac_list(k)%aperture0 = aperture(k)
     frac_list(k)%aperture  = aperture(k)
     frac_list(k)%pressure  = pressure(k)
  end do

  ! ================================================================
  ! 3. Time loop ---------------------------------------------------
  ! ================================================================
  do step = 1, nstep

     ! ----------------- Picard outer loop -------------------------
     real*8 , allocatable :: p_iter(:), p_new(:)
     if (.not.allocated(p_iter)) allocate(p_iter(NELE), p_new(NELE))
     p_iter = pressure
     do pic = 1, maxPic
        ! (a) mechanics update using current pressure guess
        do k = 1, NELE
           frac_list(k)%pressure = p_iter(k)
        end do
        call mech_update(mu_fluid, dt, ds)
        do k = 1, NELE
           aperture(k) = frac_list(k)%aperture
           if (status(k)==1) aperture(k) = max(aperture(k), w0)
        end do

        ! (b) Carter leak update for NF faces
        call update_leak_nf(C_L, step*dt, dt)

        ! (c) fluid solve over entire network (includes NF leak):
        !     we call with HF range 1:NELE for convenience; leak_hf zeros.
        real*8, allocatable :: leak_hf(:)
        if (.not.allocated(leak_hf)) allocate(leak_hf(NELE))
        leak_hf = 0.0d0
        call solve_hf_pressure(NELE, aperture, mu_fluid, dt, ds, p_iter, &
                               0, 0.0d0, leak_hf, p_new)

        if (maxval(abs(p_new - p_iter)) < tolPic*maxval(abs(p_iter)+1.d-12)) exit
        p_iter = p_new
     end do  ! Picard
     pressure = p_new

     ! ----------------- Fracture criterion ------------------------
     do ic=1,nClust
        call get_KI_tip_ddm(iEnd(ic), KI(iEnd(ic)))
     end do
     call update_fracture_lefm(NELE, KI, KIc, status)
     prev_status = status

     ! ----------------- Output ------------------------------------
     call write_state('state.json', step, pressure, aperture)
  end do

  write(*,*) 'coupled_driver finished.'
end program coupled_driver 
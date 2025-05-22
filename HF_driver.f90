program HF_driver
  !------------------------------------------------------------------
  ! Dynamic hydraulic-fracture driver (step-1 implementation)
  !   – couples the existing Poiseuille flow module with the simple
  !     linear-normal-stiffness mechanical response.
  !   – fracture propagation still uses the original tensile check.
  !   – JSON I/O and LEFM K_I will be integrated in later steps.
  !------------------------------------------------------------------
  use flow_network     ! bridge Poiseuille solver (future network)
  use mech_bridge      ! aperture from DDM compliance + NF bridge
  use fracture_crit    ! failure criterion (tensile until upgraded)
  use json_io          ! minimal JSON dumping
  use json_io_full,   only: read_config, write_state
  use fracture_utils   ! centroid loader & spacing calc (renamed)
  use ki_ddm           ! DDM-based KI computation
  use fracture_types   ! global fracture list
  implicit none

  !--------------- simulation controls ------------------------------
  integer :: NELE              ! number of quadrilateral faces along fracture mid-line
  integer :: nstep, step
  real*8  :: mu_fluid, dt, inlet_val, t0, KIc, C_L
  integer :: inlet_type
  real*8  :: ds
  real*8, allocatable :: gx(:), gy(:), gz(:)
  integer :: ierr
  integer :: k, idx_init

  !--------------- field arrays -------------------------------------
  real*8,  allocatable :: aperture(:)   ! fracture opening [m]
  real*8,  allocatable :: pressure(:)   ! fluid pressure [Pa]
  real*8,  allocatable :: sigma_n(:)    ! normal stress on each face [Pa]
  integer, allocatable :: status(:)     ! 0 = intact ; 1 = failed
  integer, allocatable :: prev_status(:)
  real*8, allocatable :: KI(:)
  real*8, allocatable :: t_open(:), leak(:)

  !--------------- local vars ---------------------------------------
  character(len=256) :: cfg_path
  logical :: use_json

  !------------------------------------------------------------------
  ! 1.  Read text-based configuration (JSON will come later)         
  !------------------------------------------------------------------
  inquire(file='config.json',exist=use_json)
  if (use_json) then
     call read_config('config.json', NELE, mu_fluid, dt, nstep, inlet_type, inlet_val, t0, KIc, C_L)
  else
     cfg_path = 'input_hf.txt'
     open(unit=10,file=cfg_path,status='old',action='read',iostat=step)
     if (step /= 0) then
        write(*,*) 'ERROR: cannot open ', trim(cfg_path)
        stop
     end if
     read(10,*) NELE, mu_fluid, dt, nstep, inlet_type, inlet_val, t0, KIc, C_L
     close(10)
  end if

  !------------------------------------------------------------------
  ! Try to compute average element spacing from output_part1.txt
  !------------------------------------------------------------------
  allocate(gx(NELE), gy(NELE), gz(NELE))
  call load_centers('../part1/output_part1.txt', NELE, gx, gy, gz, ierr)
  if (ierr == 0) then
     ds = average_spacing(NELE, gx, gy, gz)
  else
     ds = 1.0d0
  end if
  deallocate(gx, gy, gz)

  !------------------------------------------------------------------
  ! 2.  Allocate & initialise state                                   
  !------------------------------------------------------------------
  allocate(aperture(NELE), pressure(NELE), sigma_n(NELE), status(NELE))
  allocate(prev_status(NELE))
  allocate(KI(NELE))
  allocate(t_open(NELE)); t_open=-1.0d0
  allocate(leak(NELE)); leak=0.0d0
  aperture = 1.0d-3      ! 1 mm initial aperture
  pressure = 0.0d0       ! hydrostatic
  sigma_n  = 0.0d0
  status   = 0
  prev_status = 0

  !--------------------------------------------------------------
  ! Initialise global frac_list with the current HF elements
  !--------------------------------------------------------------
  if (.not.allocated(frac_list)) allocate(frac_list(NELE))
  do idx_init = 1, NELE
     frac_list(idx_init)%is_hf     = .true.
     frac_list(idx_init)%is_active = .true.
     frac_list(idx_init)%aperture0 = aperture(idx_init)
     frac_list(idx_init)%aperture  = aperture(idx_init)
     frac_list(idx_init)%pressure  = pressure(idx_init)
  end do

  !------------------------------------------------------------------
  ! 3.  Time-stepping loop                                            
  !------------------------------------------------------------------
  do step = 1, nstep

     !––– (a) mechanical update via mech_bridge –––––––––––––––––––
     do k = 1, NELE
        frac_list(k)%pressure = pressure(k)
     end do
     call mech_update(mu_fluid, dt, ds)
     do k = 1, NELE
        aperture(k) = frac_list(k)%aperture
     end do
     sigma_n = -pressure

     ! If an element is broken, ensure it has at least a base aperture
     real*8, parameter :: w0 = 1.0d-3   ! 1 mm base opening for fully broken faces
     do k = 1, NELE
        if (status(k) == 1) aperture(k) = max(aperture(k), w0)
     end do

     !––– (b) solve 1-D Poiseuille flow ––––––––––––––––––––––––––––
     ! Carter leak-off calculation
     real*8 :: current_time
     current_time = step*dt
     do k=1,NELE
        if (status(k)==1 .and. t_open(k)<0) t_open(k)=current_time-dt
        if (t_open(k)>0) then
           leak(k)=C_L/sqrt(max(current_time-t_open(k),dt))
        else
           leak(k)=0.0d0
        end if
     end do

     ! If rate-controlled inlet, iterate on inlet pressure to match rate
     if (inlet_type == 1) then
        real*8 :: Pin, Qcalc, k12
        integer :: it
        Pin = inlet_val/ (1.d-12)  ! crude initial guess
        do it=1,10
           call solve_hf_pressure(NELE, aperture, mu_fluid, dt, ds, pressure, &
                                  0, Pin, leak, pressure)
           k12 = ((aperture(1)+aperture(2))/2.d0)**3 /(12.d0*mu_fluid)
           Qcalc = k12 * (Pin - pressure(2)) / ds
           if (abs(Qcalc - inlet_val) < 1.d-6*inlet_val + 1.d-12) exit
           Pin = Pin + (inlet_val - Qcalc)*ds/k12
        end do
     else
        ! Picard fixed-point iteration between mechanics and fluid
        integer :: pic
        integer, parameter :: maxPic = 10
        real*8, parameter :: tolPic = 1.0d-4
        real*8, allocatable :: p_iter(:), p_new(:)
        if (.not.allocated(p_iter)) allocate(p_iter(NELE), p_new(NELE))
        p_iter = pressure  ! start from last step pressure

        do pic = 1, maxPic
           ! (i) mechanics with current pressure guess
           do k = 1, NELE
              frac_list(k)%pressure = p_iter(k)
           end do
           call mech_update(mu_fluid, dt, ds)
           do k = 1, NELE
              aperture(k) = frac_list(k)%aperture
              if (status(k) == 1) aperture(k) = max(aperture(k), w0)
           end do
           sigma_n = -p_iter

           ! (ii) fluid solve to update pressure
           call solve_hf_pressure(NELE, aperture, mu_fluid, dt, ds, p_iter, &
                                  inlet_type, inlet_val, leak, p_new)

           if (maxval(abs(p_new - p_iter)) < tolPic*maxval(abs(p_iter)+1.d-12)) then
              exit  ! converged
           end if
           p_iter = p_new
        end do

        pressure = p_new
     end if

     !––– (c) fracture-status update –––––––––––––––––––––––––––––––
     call get_KI_tip_ddm(NELE, KI(NELE))
     call update_fracture_lefm(NELE, KI, KIc, status)

     ! mark new broken faces traction-free
     do k = 1, NELE
        if (status(k)==1 .and. prev_status(k)==0) then
           call set_face_broken(k)
        end if
     end do
     prev_status = status

     !––– (e) simple mesh extension: add one element if the tip
     !     element has failed during this step
     if (status(NELE) == 1) then
        integer :: oldN
        oldN = NELE
        NELE  = NELE + 1

        ! Resize allocatable arrays keeping existing data
        call extend_real(aperture, oldN)
        call extend_real(pressure, oldN)
        call extend_real(sigma_n,  oldN)
        call extend_int (status,   oldN)
        call extend_real(KI,       oldN)
        call extend_real(t_open,   oldN)
        call extend_real(leak,     oldN)

        ! Initialise new face with last-known values
        aperture(NELE) = aperture(oldN)
        pressure(NELE) = pressure(oldN)
        sigma_n(NELE)  = sigma_n(oldN)
        status(NELE)   = 0
        KI(NELE)       = 0.0d0
        t_open(NELE)   = -1.0d0
        leak(NELE)     = 0.0d0

        ! Extend frac_list with a new HF element
        if (allocated(frac_list)) then
           type(frac_elem), allocatable :: tmpFL(:)
           call move_alloc(frac_list, tmpFL)
           allocate(frac_list(oldN+1))
           frac_list(1:oldN) = tmpFL
        else
           allocate(frac_list(NELE))
        end if
        frac_list(NELE)%is_hf     = .true.
        frac_list(NELE)%is_active = .true.
        frac_list(NELE)%aperture0 = aperture(NELE)
        frac_list(NELE)%aperture  = aperture(NELE)
        frac_list(NELE)%pressure  = pressure(NELE)

        ! inform ki_ddm to create geometry for the new element
        call extend_mesh(ds)
     end if

     !––– (d) write state to lightweight JSON file ––––––––––––––––
     if (use_json) then
        call write_state('state.json', step, pressure, aperture)
     else
        call dump_state(step, pressure, aperture)
     end if
  end do

  !------------------------------------------------------------------
  ! 4.  Simple output (pressure along faces)                          
  !------------------------------------------------------------------
  open(unit=20,file='pressure_out.txt',status='replace',action='write')
  do step = 1, NELE
     write(20,'(I8,1x,ES14.5)') step, pressure(step)
  end do
  close(20)

  write(*,*) 'HF_driver finished. Results in pressure_out.txt'
end program HF_driver

contains
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
end program HF_driver 
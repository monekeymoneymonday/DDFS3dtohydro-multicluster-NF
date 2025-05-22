!*********************************************************************
! Example driver to couple mechanical solver with hydraulic fracture
!*********************************************************************
program HF_main
  use flow_hf
  use fracture_crit
  use mechanics_hf
  use subs_fs3d_qua_num
  implicit none

  integer :: NELE, step, nstep
  real*8  :: mu_fluid, dt, inj_rate, t0, krock
  real*8, allocatable :: aperture(:), pressure(:), sigma_n(:)
  integer, allocatable :: status(:)

  !-----------------------------------------------------------------
  ! Input section
  open(20,file='input_hf.txt')
  read(20,*) NELE, mu_fluid, dt, nstep, inj_rate, t0, krock
  close(20)

  allocate(aperture(NELE), pressure(NELE), sigma_n(NELE), status(NELE))
  aperture = 1.0e-3      ! initial fracture aperture [m]
  pressure = 0.0         ! pore-pressure field [Pa]
  sigma_n  = 0.0         ! normal stress on each element [Pa]
  status   = 0           ! 0 = intact, 1 = broken

  do step = 1, nstep
     !––– Mechanics solve (updates aperture & σₙ) ––––––––––––––––––––
     call update_mechanics(NELE, pressure, krock, aperture, sigma_n)

     !––– Fluid-flow solve with updated aperture –––––––––––––––––––––
     call solve_hf_pressure(NELE, aperture, mu_fluid, dt, pressure,     &
                            inj_rate, pressure)

     !––– Fracture-propagation criterion –––––––––––––––––––––––––––––
     call update_fracture_status(NELE, sigma_n, pressure, t0, status)
  end do

  !-----------------------------------------------------------------
  ! Output
  open(30,file='pressure_out.txt')
  do step = 1, NELE
     write(30,'(I8,ES12.4)') step, pressure(step)
  end do
  close(30)
end program HF_main

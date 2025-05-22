module flow_hf
  implicit none
contains
  !---------------------------------------------------------------
  ! Thomas algorithm for a tridiagonal linear system  (A x = rhs)
  !---------------------------------------------------------------
  subroutine solve_tridiag(n, a, b, c, rhs, x)
    integer, intent(in)    :: n
    real*8, intent(inout)  :: a(n), b(n), c(n), rhs(n)
    real*8, intent(out)    :: x(n)
    integer :: i

    do i = 2, n
       a(i)  = a(i) / b(i-1)
       b(i)  = b(i) - a(i) * c(i-1)
       rhs(i)= rhs(i) - a(i) * rhs(i-1)
    end do

    x(n) = rhs(n) / b(n)
    do i = n-1, 1, -1
       x(i) = (rhs(i) - c(i) * x(i+1)) / b(i)
    end do
  end subroutine solve_tridiag

  !---------------------------------------------------------------
  ! Backward-Euler Poiseuille flow solver along a 1-D fracture
  !---------------------------------------------------------------
  subroutine solve_hf_pressure(num_faces, aperture, mu, dt, ds, p_old, inlet_type, inlet_val, leak, p_new)
    integer, intent(in)  :: num_faces
    real*8 , intent(in)  :: aperture(num_faces)
    real*8 , intent(in)  :: mu, dt
    real*8 , intent(in), optional :: ds        ! ds = face-to-face spacing [m]
    real*8 , intent(in)  :: p_old(num_faces)
    integer, intent(in)  :: inlet_type         ! 0=Dirichlet P, 1=rate Q
    real*8 , intent(in)  :: inlet_val          ! pressure [Pa] or rate [m3/s]
    real*8 , intent(in)  :: leak(num_faces)    ! leak-off volume rate [m^2/s]
    real*8 , intent(out) :: p_new(num_faces)

    real*8 :: a(num_faces), b(num_faces), c(num_faces), rhs(num_faces)
    real*8 :: kL, kR, ds_loc
    integer :: i

    ! Determine spacing
    if (present(ds)) then
       ds_loc = ds
    else
       ds_loc = 1.0d0
    end if

    ! Interior nodes
    do i = 2, num_faces-1
       kL      = ((aperture(i-1) + aperture(i  )) / 2.d0)**3 / (12.d0 * mu)
       kR      = ((aperture(i  ) + aperture(i+1)) / 2.d0)**3 / (12.d0 * mu)
       a(i)    = -dt * kL / (ds_loc * ds_loc)
       c(i)    = -dt * kR / (ds_loc * ds_loc)
       b(i)    = 1.d0 - a(i) - c(i)
       rhs(i)  = p_old(i) - dt*leak(i)/max(aperture(i),1d-9)
    end do

    ! Inlet boundary treatment
    if (inlet_type == 0) then  ! Dirichlet pressure
       a(1)=0.d0; c(1)=0.d0; b(1)=1.d0
       rhs(1)=inlet_val
    else                        ! Rate-controlled source (simplified)
       a(1)=0.d0; c(1)=0.d0; b(1)=1.d0
       rhs(1)=p_old(1) + dt*inlet_val/(max(aperture(1),1d-9)*ds_loc)
    end if

    ! Fracture tip (zero-flux Neumann)
    kL            = ((aperture(num_faces-1) + aperture(num_faces)) / 2.d0)**3 &
                    / (12.d0 * mu)
    a(num_faces)  = -dt * kL / (ds_loc * ds_loc)
    b(num_faces)  = 1.d0 - a(num_faces)
    c(num_faces)  = 0.d0
    rhs(num_faces)= p_old(num_faces) - dt*leak(num_faces)/max(aperture(num_faces),1d-9)

    ! Solve the tridiagonal system
    call solve_tridiag(num_faces, a, b, c, rhs, p_new)
  end subroutine solve_hf_pressure
end module flow_hf

module flow_hf_multi
  !--------------------------------------------------------------------
  ! Thin wrapper that solves the 1-D Poiseuille flow problem for
  ! several fracture clusters in one call. Each cluster is assumed
  ! to be represented by a *contiguous* range of quadrilateral faces.
  ! The heavy lifting is still performed by the original
  ! flow_hf::solve_hf_pressure routine.
  !--------------------------------------------------------------------
  use flow_hf, only : solve_hf_pressure
  implicit none
  private
  public :: solve_multi_pressure
contains
  !------------------------------------------------------------------
  ! Loop over clusters and call the existing solver on each range.
  !------------------------------------------------------------------
  subroutine solve_multi_pressure(nClust, iStart, iEnd,                 &
                                  aperture, mu, dt, ds,                &
                                  inletType, inletVal, leak, p)
    integer, intent(in)  :: nClust                ! number of clusters
    integer, intent(in)  :: iStart(nClust)        ! first element index of each cluster (1-based)
    integer, intent(in)  :: iEnd(nClust)          ! last  element index of each cluster (1-based)

    real*8 , intent(in)  :: aperture(:)           ! total-field aperture array
    real*8 , intent(in)  :: mu                    ! fluid viscosity [Pa·s]
    real*8 , intent(in)  :: dt                    ! time-step [s]
    real*8 , intent(in)  :: ds                    ! element spacing [m]
    integer, intent(in)  :: inletType(nClust)     ! 0 = Dirichlet P ; 1 = rate Q
    real*8 , intent(in)  :: inletVal(nClust)      ! Pa or m³/s (same size as clusters)
    real*8 , intent(in)  :: leak(:)               ! Carter leak-off for each element [m²/s]
    real*8 , intent(inout) :: p(:)                ! pressure field (in/out) [Pa]

    integer :: ic, nLoc, i1, i2

    do ic = 1, nClust
       i1   = iStart(ic)
       i2   = iEnd(ic)
       nLoc = i2 - i1 + 1
       if (nLoc <= 0) cycle   ! safety guard

       call solve_hf_pressure(nLoc,                         &
                              aperture(i1:i2),              &
                              mu, dt, ds,                   &
                              p(i1:i2),                     &
                              inletType(ic), inletVal(ic),  &
                              leak(i1:i2),                  &
                              p(i1:i2))
    end do
  end subroutine solve_multi_pressure

end module flow_hf_multi 
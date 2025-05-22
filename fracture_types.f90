module fracture_types
  implicit none
  private

  ! Double precision kind parameter (15 decimal digits)
  integer, parameter, public :: dp = selected_real_kind(15, 307)

  !------------------------------------------------------------------
  ! Unified fracture element – represents both hydraulic and natural
  ! fractures.  Geometry is stored as the 4 corner points of a
  ! quadrilateral displacement–discontinuity element (same layout as
  ! existing FSM3D-Q formulation).
  !------------------------------------------------------------------
  type, public :: frac_elem
     ! Geometry (global coordinates)
     real(dp) :: x(4) = 0.0_dp, y(4) = 0.0_dp, z(4) = 0.0_dp
     real(dp) :: normal(3)   = 0.0_dp  ! unit normal vector

     ! Mechanical parameters
     real(dp) :: Kn = 0.0_dp      ! normal stiffness [Pa/m]
     real(dp) :: Ks = 0.0_dp      ! shear  stiffness [Pa/m]
     real(dp) :: aperture0 = 0.0_dp  ! initial hydraulic aperture [m]
     real(dp) :: mu = 0.0_dp         ! friction coefficient μ
     real(dp) :: cohesion = 0.0_dp   ! cohesion [Pa]
     real(dp) :: tensile_strength = 0.0_dp ! tensile strength [Pa]

     ! Dynamic state (updated each simulation step)
     logical  :: is_active = .false.   ! true → element currently conductive
     logical  :: is_hf     = .false.   ! true → element belongs to HF
     logical  :: can_grow  = .false.   ! true → this is a propagating tip
     real(dp) :: aperture  = 0.0_dp    ! current aperture [m]
     real(dp) :: pressure  = 0.0_dp    ! current fluid pressure [Pa]
     real(dp) :: shear_disp = 0.0_dp   ! accumulated shear displacement [m]

     ! Leak-off tracking (Carter)
     real(dp) :: t_open   = -1.0_dp   ! time when element became active [s]; <0 ⇒ not yet
     real(dp) :: leak     = 0.0_dp    ! current Carter leak-off rate [m^2/s]

     ! Connectivity (to be filled when meshing)
     integer, allocatable :: neighbors(:)  ! indices of adjacent elements
  end type frac_elem

  !------------------------------------------------------------------
  ! Global container holding *all* fracture elements (HF + NF).
  ! Subsequent modules (mechanics, flow, etc.) will operate on this
  ! list so that both fracture types are treated in a unified manner.
  !------------------------------------------------------------------
  type(frac_elem), allocatable, public :: frac_list(:)

end module fracture_types 
module nf_types
  !------------------------------------------------------------------
  ! Natural-fracture data catalogue and lightweight text reader.
  ! Format of nf_catalog.txt (one fracture per line):
  !   id   x1 y1 z1  x2 y2 z2  x3 y3 z3  x4 y4 z4   aperture0  Kn  Ks  Cohesion  mu  Tensile
  ! Units: coordinates in metres, aperture [m], stiffness [Pa/m], cohesion & tensile strength [Pa].
  ! Lines beginning with # are ignored.
  !------------------------------------------------------------------
  use fracture_types, only: dp
  implicit none
  private
  public :: NF, nf_list, read_nf_catalog

  type :: NF
     integer :: id = 0
     real(dp) :: x(4) = 0.0_dp, y(4)=0.0_dp, z(4)=0.0_dp
     real(dp) :: aperture0 = 1.0e-4_dp
     real(dp) :: Kn = 1.0e12_dp          ! normal stiffness
     real(dp) :: Ks = 1.0e12_dp          ! shear  stiffness
     real(dp) :: cohesion = 0.0_dp       ! cohesion [Pa]
     real(dp) :: mu = 0.6_dp             ! friction coefficient
     real(dp) :: tensile = 0.0_dp        ! tensile strength [Pa]
     integer  :: state = 0               ! 0=sealed, 1=open, 2=slipped
  end type NF

  type(NF), allocatable, public :: nf_list(:)

contains

  !----------------------------------------------------------------
  subroutine read_nf_catalog(filename)
    character(len=*), intent(in) :: filename
    integer :: ios, nLines, i
    character(len=1024) :: line

    ! First pass: count fractures
    nLines = 0
    open(unit=99,file=filename,status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'read_nf_catalog: cannot open catalog'
    do
       read(99,'(A)',iostat=ios) line
       if (ios/=0) exit
       if (trim(line)=='' .or. line(1:1)=='#') cycle
       nLines = nLines + 1
    end do
    close(99)
    if (nLines <= 0) return

    allocate(nf_list(nLines))

    ! Second pass: parse values
    open(unit=99,file=filename,status='old',action='read')
    i = 0
    do
       read(99,'(A)',iostat=ios) line
       if (ios/=0) exit
       if (trim(line)=='' .or. line(1:1)=='#') cycle
       i = i + 1
       read(line,*) nf_list(i)%id, &
            nf_list(i)%x(1), nf_list(i)%y(1), nf_list(i)%z(1), &
            nf_list(i)%x(2), nf_list(i)%y(2), nf_list(i)%z(2), &
            nf_list(i)%x(3), nf_list(i)%y(3), nf_list(i)%z(3), &
            nf_list(i)%x(4), nf_list(i)%y(4), nf_list(i)%z(4), &
            nf_list(i)%aperture0, nf_list(i)%Kn, nf_list(i)%Ks, &
            nf_list(i)%cohesion, nf_list(i)%mu, nf_list(i)%tensile
    end do
    close(99)
  end subroutine read_nf_catalog

end module nf_types 
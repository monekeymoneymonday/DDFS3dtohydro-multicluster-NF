module ki_ddm
  use subs_fs3d_qua_num
  implicit none
  logical :: loaded = .false.
  integer :: nEle, nNod
  real*8  :: E_mat, nu_mat
  integer, allocatable :: ELE(:,:)
  real*8 , allocatable :: NOD(:,:)
  real*8 , allocatable :: Pvec(:,:)
  logical, public :: Cnn_dirty = .true.   ! signal compliance needs rebuild
contains
  !----------------------------------------------------------------------
  subroutine init_ki_data()
    character(len=*), parameter :: finp = '../part1/input_part1.txt'
    character(len=*), parameter :: fp   = '../part1/P.txt'
    integer :: ios, i, dummyI
    real*8 :: dummyR
    character(len=256) :: title

    open(unit=21,file=finp,status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'init_ki_data: cannot open input_part1.txt'

    read(21,'(A)') title
    read(21,*) E_mat, nu_mat, dummyR, dummyR  ! KS, KN ignored
    read(21,*) dummyR, dummyR, dummyR, dummyR, dummyR, dummyR   ! pxx..pyz
    read(21,*) nEle, nNod

    allocate(ELE(nEle,4), NOD(nNod,3), Pvec(nEle,3))

    do i=1, nEle
       read(21,*) dummyI, ELE(i,1:4)
    end do
    do i=1, nNod
       read(21,*) dummyI, NOD(i,1:3)
    end do
    ! skip boundary section
    ! we don't need SOB/VOB for KI
    close(21)

    ! read fictitious stresses
    open(unit=22,file=fp,status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'init_ki_data: cannot open P.txt'
    do i=1, nEle
       read(22,*) dummyI, Pvec(i,1:3)
    end do
    close(22)

    loaded = .true.
  end subroutine init_ki_data

  !----------------------------------------------------------------------
  subroutine get_KI_tip_ddm(elem_id, KI)
    integer, intent(in)  :: elem_id
    real*8 , intent(out) :: KI

    real*8, parameter :: pi = 3.141592653589793d0
    real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    real*8 :: IEV(3,3), COEFF(3,3)
    real*8 :: sigma_z
    integer :: M

    if (.not.loaded) call init_ki_data()

    if (elem_id < 1 .or. elem_id > nEle) then
       KI = 0.0d0
       return
    end if

    ! node coords
    x1 = NOD( ELE(elem_id,1), 1 ); y1 = NOD( ELE(elem_id,1), 2 ); z1 = NOD( ELE(elem_id,1), 3 )
    x2 = NOD( ELE(elem_id,2), 1 ); y2 = NOD( ELE(elem_id,2), 2 ); z2 = NOD( ELE(elem_id,2), 3 )
    x3 = NOD( ELE(elem_id,3), 1 ); y3 = NOD( ELE(elem_id,3), 2 ); z3 = NOD( ELE(elem_id,3), 3 )
    x4 = NOD( ELE(elem_id,4), 1 ); y4 = NOD( ELE(elem_id,4), 2 ); z4 = NOD( ELE(elem_id,4), 3 )

    ! local basis
    call convert(x1,y1,z1,x2,y2,z2,x3,y3,z3, IEV)

    M = 18
    call stress_singular_fs3d_qua_num(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                                      E_mat, nu_mat, M, IEV, COEFF)

    ! stress ahead of tip (local z-direction normal): combine with fictitious stresses
    sigma_z = COEFF(3,1)*Pvec(elem_id,1) + COEFF(3,2)*Pvec(elem_id,2) + COEFF(3,3)*Pvec(elem_id,3)

    ! K_I = sigma * sqrt(pi * a_eff) ; choose a_eff ~ element length / 2
    real*8 :: a_eff, dx,dy,dz
    dx = 0.5d0*((x2 - x1)+(x3 - x4))
    dy = 0.5d0*((y2 - y1)+(y3 - y4))
    dz = 0.5d0*((z2 - z1)+(z3 - z4))
    a_eff = 0.5d0*sqrt(dx*dx + dy*dy + dz*dz)

    KI = sigma_z * sqrt(pi * max(a_eff,1.0d-6))
  end subroutine get_KI_tip_ddm

  !----------------------------------------------------------------------
  ! Extend mesh by cloning the current tip element and translating it
  ! a distance ds along the local x-axis (propagation direction).
  !----------------------------------------------------------------------
  subroutine extend_mesh(ds)
    real*8, intent(in) :: ds

    integer :: oldEle, i, oldNod
    real*8 :: IEV(3,3), shift(3)
    integer, allocatable :: tmpI(:,:)
    real*8 , allocatable :: tmpR(:,:), tmpP(:,:)

    if (.not.loaded) call init_ki_data()

    oldEle = nEle
    oldNod = nNod

    ! compute local basis of tip element
    call convert( NOD(ELE(oldEle,1),1), NOD(ELE(oldEle,1),2), NOD(ELE(oldEle,1),3), &
                 NOD(ELE(oldEle,2),1), NOD(ELE(oldEle,2),2), NOD(ELE(oldEle,2),3), &
                 NOD(ELE(oldEle,3),1), NOD(ELE(oldEle,3),2), NOD(ELE(oldEle,3),3), IEV )

    shift = ds * IEV(1,:)  ! global propagation vector

    ! append 4 new nodes (duplicates shifted)
    allocate(tmpR(nNod,3)); tmpR = NOD; call move_alloc(tmpR, NOD)
    allocate(tmpR(4,3)); tmpR = 0.0d0  ! new nodes container
    do i=1,4
       tmpR(i,1) = NOD( ELE(oldEle,i), 1 ) + shift(1)
       tmpR(i,2) = NOD( ELE(oldEle,i), 2 ) + shift(2)
       tmpR(i,3) = NOD( ELE(oldEle,i), 3 ) + shift(3)
    end do
    ! increase NOD size
    allocate(tmpR(nNod+4,3))
    tmpR(1:nNod, :) = NOD
    tmpR(nNod+1:nNod+4, :) = reshape(tmpR, [4,3])
    call move_alloc(tmpR, NOD)
    nNod = nNod + 4

    ! append new element connectivity
    allocate(tmpI(nEle,4)); tmpI = ELE; call move_alloc(tmpI, ELE)
    allocate(tmpI(nEle+1,4))
    tmpI(1:nEle, :) = ELE
    do i=1,4
       tmpI(nEle+1,i) = oldNod + i
    end do
    call move_alloc(tmpI, ELE)
    nEle = nEle + 1

    ! append fictitious stress zeros for new element
    allocate(tmpP(size(Pvec,1),3)); tmpP = Pvec; call move_alloc(tmpP, Pvec)
    allocate(tmpP(nEle,3))
    tmpP(1:nEle-1,:) = Pvec
    tmpP(nEle,:) = 0.0d0
    call move_alloc(tmpP, Pvec)
  end subroutine extend_mesh

  !----------------------------------------------------------------------
  subroutine set_face_broken(id)
    integer, intent(in) :: id
    if (.not.loaded) call init_ki_data()
    if (id<1 .or. id>nEle) return
    Pvec(id,1:3) = 0.0d0
    Cnn_dirty = .true.   ! signal compliance needs rebuild
  end subroutine set_face_broken

  !----------------------------------------------------------------------
  ! Placeholder: insert a new element *after* a given element index.
  ! A full implementation would shift internal arrays so that the new
  ! element appears at position id+1, preserving cluster contiguity.
  ! For now we simply fall back to extend_mesh (appending at the end),
  ! which means the ordering may no longer be contiguous.  The driving
  ! code must account for this or the routine should be upgraded later.
  !----------------------------------------------------------------------
  subroutine extend_mesh_insert(id, ds)
    integer, intent(in) :: id
    real*8 , intent(in) :: ds

    integer :: i, newNodes(4)
    real*8 :: IEV(3,3), shift(3)
    integer, allocatable :: tmpI(:,:)
    real*8 , allocatable :: tmpR(:,:)

    if (.not.loaded) call init_ki_data()

    if (id < 1 .or. id > nEle) then
       call extend_mesh(ds)   ! fallback
       return
    end if

    if (id == nEle) then
       call extend_mesh(ds)
       return
    end if

    !-----------------------------------------------------------
    ! (1) create 4 new nodes shifted from element 'id'
    !-----------------------------------------------------------
    call convert( NOD(ELE(id,1),1), NOD(ELE(id,1),2), NOD(ELE(id,1),3), &
                 NOD(ELE(id,2),1), NOD(ELE(id,2),2), NOD(ELE(id,2),3), &
                 NOD(ELE(id,3),1), NOD(ELE(id,3),2), NOD(ELE(id,3),3), IEV )

    shift = ds * IEV(1,:)

    allocate(tmpR(nNod,3)); tmpR = NOD; call move_alloc(tmpR, NOD)
    allocate(tmpR(4,3))
    do i=1,4
       tmpR(i,1) = NOD( ELE(id,i),1 ) + shift(1)
       tmpR(i,2) = NOD( ELE(id,i),2 ) + shift(2)
       tmpR(i,3) = NOD( ELE(id,i),3 ) + shift(3)
       newNodes(i) = nNod + i
    end do
    allocate(tmpR(nNod+4,3))
    tmpR(1:nNod,:) = NOD
    tmpR(nNod+1:nNod+4,:) = reshape(tmpR, [4,3])
    call move_alloc(tmpR, NOD)
    nNod = nNod + 4

    !-----------------------------------------------------------
    ! (2) enlarge ELE & Pvec with insertion after 'id'
    !-----------------------------------------------------------
    allocate(tmpI(nEle+1,4))
    tmpI(1:id,:) = ELE(1:id,:)
    tmpI(id+1,:) = newNodes
    tmpI(id+2:nEle+1,:) = ELE(id+1:nEle,:)
    call move_alloc(tmpI, ELE)

    allocate(tmpR(nEle+1,3))
    tmpR(1:id,:) = Pvec(1:id,:)
    tmpR(id+1,:) = 0.0d0
    tmpR(id+2:nEle+1,:) = Pvec(id+1:nEle,:)
    call move_alloc(tmpR, Pvec)

    nEle = nEle + 1
  end subroutine extend_mesh_insert
end module ki_ddm 
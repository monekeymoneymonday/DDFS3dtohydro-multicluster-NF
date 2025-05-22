module ddm_mechanics
  ! Computes fracture apertures from pressure using the displacement-
  ! discontinuity compliance matrix that already exists in the codebase.
  use subs_fs3d_qua_num
  use ki_ddm, only: ELE, NOD, nEle, loaded, init_ki_data, Cnn_dirty, convert, E_mat, nu_mat, Pvec
  implicit none
  private

  real*8, allocatable :: Cnn(:,:)  ! normal compliance matrix
  integer :: cached_ele = -1       ! number of elements that Cnn matches
  integer, parameter :: ver = 1

  public :: aperture_from_pressure, get_elem_stress

contains

  !------------------------------------------------------------------
  subroutine build_compliance()
    character(len=*), parameter :: cacheFile='Cnn.bin'
    logical :: file_exists
    integer :: ios_bin
    integer :: i, j, M
    real*8 :: IEV(3,3), JEV(3,3), COEFF(3,3)
    real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    real*8 :: xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4

    if (.not.loaded) call init_ki_data()   ! ensure geometry in ki_ddm

    ! try to load cache first
    inquire(file=cacheFile,exist=file_exists)
    if (file_exists .and. cached_ele/=nEle) then
       open(99,file=cacheFile,form='unformatted',access='stream',iostat=ios_bin)
       if (ios_bin==0) then
         read(99) i ! version
         read(99) j ! stored elements
         if (i==ver .and. j==nEle) then
            if (allocated(Cnn)) deallocate(Cnn)
            allocate(Cnn(nEle,nEle))
            read(99) Cnn
            close(99)
            cached_ele = nEle
            return
         end if
         close(99)
       end if
    end if

    if (allocated(Cnn)) deallocate(Cnn)
    allocate(Cnn(nEle,nEle))

    do i = 1, nEle
      x1=NOD(ELE(i,1),1); y1=NOD(ELE(i,1),2); z1=NOD(ELE(i,1),3)
      x2=NOD(ELE(i,2),1); y2=NOD(ELE(i,2),2); z2=NOD(ELE(i,2),3)
      x3=NOD(ELE(i,3),1); y3=NOD(ELE(i,3),2); z3=NOD(ELE(i,3),3)
      call convert(x1,y1,z1,x2,y2,z2,x3,y3,z3,IEV)
      do j = 1, nEle
        xx1=NOD(ELE(j,1),1); yy1=NOD(ELE(j,1),2); zz1=NOD(ELE(j,1),3)
        xx2=NOD(ELE(j,2),1); yy2=NOD(ELE(j,2),2); zz2=NOD(ELE(j,2),3)
        xx3=NOD(ELE(j,3),1); yy3=NOD(ELE(j,3),2); zz3=NOD(ELE(j,3),3)
        xx4=NOD(ELE(j,4),1); yy4=NOD(ELE(j,4),2); zz4=NOD(ELE(j,4),3)
        call convert(xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,JEV)
        M = 18; if (i /= j) M = 5
        if (i == j) then
          call dis_singular_fs3d_qua_num(xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4, &
                                         1d0,0.25d0,M,IEV,COEFF)
        else
          call dis_regular_fs3d_qua_num(xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3,xx4,yy4,zz4, &
                                         0d0,0d0,0d0,1d0,0.25d0,M,IEV,JEV,COEFF)
        end if
        Cnn(i,j) = COEFF(3,3)
      end do
    end do
    cached_ele = nEle

    ! save cache
    open(99,file=cacheFile,form='unformatted',access='stream',status='replace')
    write(99) ver
    write(99) nEle
    write(99) Cnn
    close(99)
  end subroutine build_compliance

  !------------------------------------------------------------------
  subroutine aperture_from_pressure(p, w)
    ! Compute aperture from pressure using compliance
    real*8, intent(in)  :: p(:)
    real*8, intent(out) :: w(:)

    if (nEle /= cached_ele .or. Cnn_dirty) then
       call build_compliance()
       Cnn_dirty = .false.
    end if

    if (.not.allocated(w)) allocate(w(nEle))
    w = matmul(Cnn, -p)
  end subroutine aperture_from_pressure

  !------------------------------------------------------------------
  ! Compute Cartesian stress tensor at the centroid of element `idx`.
  ! The stress is obtained by combining the pre-computed influence
  ! coefficients (returned by stress_singular_fs3d_qua_num) with the
  ! fictitious stress vector stored in ki_ddm:Pvec.
  !------------------------------------------------------------------
  subroutine get_elem_stress(idx, sigma_cart)
    use fracture_types , only: frac_list, dp
    implicit none

    integer, intent(in)          :: idx
    real(dp), intent(out)        :: sigma_cart(3,3)

    real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
    real*8 :: IEV(3,3), Ccoeff(3,3)
    real*8 :: X_local(3)
    integer :: col, Mpts

    ! Ensure mechanical data loaded
    if (.not.loaded) call init_ki_data()

    ! Guard against invalid index or missing fracture list
    if (.not.allocated(frac_list) .or. idx < 1 .or. idx > size(frac_list)) then
       sigma_cart = 0.0_dp
       return
    end if

    ! Geometry of the selected element (global coordinates)
    x1 = frac_list(idx)%x(1); y1 = frac_list(idx)%y(1); z1 = frac_list(idx)%z(1)
    x2 = frac_list(idx)%x(2); y2 = frac_list(idx)%y(2); z2 = frac_list(idx)%z(2)
    x3 = frac_list(idx)%x(3); y3 = frac_list(idx)%y(3); z3 = frac_list(idx)%z(3)
    x4 = frac_list(idx)%x(4); y4 = frac_list(idx)%y(4); z4 = frac_list(idx)%z(4)

    ! Local orthonormal basis of the element
    call convert(x1,y1,z1, x2,y2,z2, x3,y3,z3, IEV)

    ! Influence coefficients for stresses at centroid (singular case)
    Mpts = 18
    call stress_singular_fs3d_qua_num(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                                       E_mat, nu_mat, Mpts, IEV, Ccoeff)

    ! Fictitious stresses acting on the element (Px, Py, Pz)
    if (idx <= size(Pvec,1)) then
       X_local = Pvec(idx,1:3)
    else
       X_local = 0.0d0
    end if

    ! Combine coefficients with fictitious stresses.
    sigma_cart = 0.0_dp
    do col = 1, 3
       sigma_cart(:,col) = Ccoeff(:,col) * X_local(col)
    end do
  end subroutine get_elem_stress

end module ddm_mechanics 
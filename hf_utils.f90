module hf_utils
  implicit none
contains
  !----------------------------------------------------------------------
  ! Load centroid coordinates written by FSM3D_qua_part1 (output_part1.txt)
  ! The file has header lines followed by rows:
  !   i  Gx  Gy  Gz  Px  Py  Pz
  ! This routine extracts Gx, Gy, Gz for i = 1..num
  !----------------------------------------------------------------------
  subroutine load_centers(filename, num, gx, gy, gz, ierr)
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: num
    real*8,           intent(out) :: gx(num), gy(num), gz(num)
    integer,          intent(out) :: ierr

    character(len=256) :: line
    integer :: ios, idx, count
    real*8  :: x, y, z, dummy1, dummy2, dummy3

    ierr  = 0
    gx    = 0.0d0
    gy    = 0.0d0
    gz    = 0.0d0
    count = 0

    open(unit=98,file=filename,status='old',action='read',iostat=ios)
    if (ios /= 0) then
       ierr = 1
       return
    end if

    do
       read(98,'(A)',iostat=ios) line
       if (ios /= 0) exit   ! EOF

       ! try to parse numeric line
       read(line,*,iostat=ios) idx, x, y, z, dummy1, dummy2, dummy3
       if (ios == 0) then
          if (idx >= 1 .and. idx <= num) then
             gx(idx) = x
             gy(idx) = y
             gz(idx) = z
             count   = count + 1
          end if
       end if
    end do
    close(98)

    if (count /= num) ierr = 2
  end subroutine load_centers

  !----------------------------------------------------------------------
  ! Compute the average spacing between consecutive centroids
  !----------------------------------------------------------------------
  function average_spacing(num, gx, gy, gz) result(ds)
    integer, intent(in) :: num
    real*8 , intent(in) :: gx(num), gy(num), gz(num)
    real*8 :: ds, dist
    integer :: i

    if (num <= 1) then
       ds = 1.0d0
       return
    end if

    ds = 0.0d0
    do i = 1, num-1
       dist = sqrt( (gx(i+1)-gx(i))**2 + (gy(i+1)-gy(i))**2 + (gz(i+1)-gz(i))**2 )
       ds = ds + dist
    end do
    ds = ds / real(num-1)
  end function average_spacing
end module hf_utils 
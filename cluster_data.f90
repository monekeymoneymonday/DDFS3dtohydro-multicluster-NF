module cluster_data
  !--------------------------------------------------------------------
  ! Provide utilities to read simple text configuration files that
  ! specify multiple hydraulic-fracture clusters.  The format assumed
  ! (one key per line, free formatting) is for example:
  !
  !   nCluster            3
  !   NELE_per_cluster    40 35 55
  !   inletType           1 0 1
  !   inletVal            1.5e-3  5.0e6  1.0e-3
  !
  ! Keys are case-insensitive; lines starting with # are ignored.
  !--------------------------------------------------------------------
  implicit none
  private
  public :: read_cluster_config
contains
  !------------------------------------------------------------------
  subroutine read_cluster_config(file, nClust, iStart, iEnd,            &
                                 inletType, inletVal, NELE_total)
    character(len=*), intent(in)          :: file
    integer,            intent(out)       :: nClust
    integer, allocatable, intent(out)     :: iStart(:), iEnd(:)
    integer, allocatable, intent(out)     :: inletType(:)
    real*8 , allocatable, intent(out)     :: inletVal(:)
    integer,            intent(out)       :: NELE_total

    !---------------------------------------------
    ! local buffers
    character(len=256) :: line, key
    character(len=512) :: values
    integer            :: ios, i, pos
    integer, allocatable :: nElePerClust(:)

    nClust = -1
    allocate(nElePerClust(0))     ! dummy initial allocation

    open(unit=15,file=file,status='old',action='read',iostat=ios)
    if (ios /= 0) then
       write(*,*) 'ERROR: cannot open ',trim(file)
       stop
    end if

    do
       read(15,'(A)',iostat=ios) line
       if (ios /= 0) exit
       if (trim(line) == '') cycle
       if (line(1:1) == '#') cycle

       ! split into key and remainder
       pos = index(line, ' ')
       if (pos == 0) cycle
       key    = adjustl( lower(line(1:pos-1)) )
       values = adjustl( line(pos:) )

       select case (trim(key))
       case ('ncluster')
          read(values,*) nClust
       case ('nele_per_cluster')
          if (nClust <= 0) then
             write(*,*) 'read_cluster_config: nCluster must appear before NELE_per_cluster'
             stop
          end if
          if (allocated(nElePerClust)) deallocate(nElePerClust)
          allocate(nElePerClust(nClust))
          read(values,*) (nElePerClust(i), i=1,nClust)
       case ('inlettype')
          if (nClust <= 0) cycle
          if (allocated(inletType)) deallocate(inletType)
          allocate(inletType(nClust))
          read(values,*) (inletType(i), i=1,nClust)
       case ('inletval')
          if (nClust <= 0) cycle
          if (allocated(inletVal)) deallocate(inletVal)
          allocate(inletVal(nClust))
          read(values,*) (inletVal(i), i=1,nClust)
       end select
    end do
    close(15)

    ! sanity checks
    if (nClust <= 0) then
       write(*,*) 'read_cluster_config: missing or invalid nCluster'
       stop
    end if
    if (.not.allocated(nElePerClust)) then
       write(*,*) 'read_cluster_config: missing NELE_per_cluster'
       stop
    end if
    if (.not.allocated(inletType)) then
       allocate(inletType(nClust)); inletType=0
    end if
    if (.not.allocated(inletVal)) then
       allocate(inletVal(nClust)); inletVal=0.0d0
    end if

    ! compute iStart/iEnd and total number of faces
    allocate(iStart(nClust), iEnd(nClust))
    iStart(1) = 1
    do i=1,nClust
       if (i>1) iStart(i) = iEnd(i-1)+1
       iEnd(i)   = iStart(i) + nElePerClust(i) - 1
    end do
    NELE_total = iEnd(nClust)

    deallocate(nElePerClust)
  end subroutine read_cluster_config

  !------------------------------------------------------------------
  pure function lower(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str))      :: out
    integer :: k
    do k=1,len(str)
       select case(str(k:k))
       case('A':'Z')
          out(k:k)=achar(iachar(str(k:k))+32)
       case default
          out(k:k)=str(k:k)
       end select
    end do
  end function lower
end module cluster_data 
module json_io
  !----------------------------------------------------------------------
  ! Minimal JSON writer (plain text) – placeholder until we link to
  ! JSON-Fortran.  Generates a simple file  state_<step>.json  that can
  ! be parsed by any JSON reader.
  !----------------------------------------------------------------------
  implicit none
contains
  subroutine dump_state(step, pressure, aperture)
    integer, intent(in) :: step
    real*8 , intent(in) :: pressure(:), aperture(:)
    integer :: i, n
    character(len=64) :: fname

    n = size(pressure)
    write(fname,'(A,I0.5,A)') 'state_', step, '.json'

    open(unit=91,file=fname,status='replace',action='write')
    write(91,'(A)') '{'
    write(91,'(A,I0,A)') '  "step": ', step, ','
    write(91,'(A)') '  "pressure": ['
    do i = 1, n-1
       write(91,'(ES16.8,A)') pressure(i), ','
    end do
    write(91,'(ES16.8)') pressure(n)
    write(91,'(A)') '  ],'
    write(91,'(A)') '  "aperture": ['
    do i = 1, n-1
       write(91,'(ES16.8,A)') aperture(i), ','
    end do
    write(91,'(ES16.8)') aperture(n)
    write(91,'(A)') '  ]'
    write(91,'(A)') '}'
    close(91)
  end subroutine dump_state
end module json_io 
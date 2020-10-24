!=======================================================================
module my_pio_routines

  use piolibkind
  use piolib
  use healpix_types
  use misc_utils, only: fatal_error

  implicit none

  private
  public :: mypioreadmap, pio_datatype, pio_objecttype

  interface mypioreadmap
     module procedure mypioreadmap_sp, mypioreadmap_dp
  end interface

  integer, private, parameter :: MYLEN = DMCPIOSTRINGMAXLEN ! AUG08-v05 and higher

contains

  !---------------------------
  function pio_objecttype(file) result (outstr)

    character(len=MYLEN) :: outstr
    character(len=*), intent(in) :: file
    integer(piolong) :: myerr

    CHARACTER(LEN=MYLEN)  :: TOITYPE, DATATYPE, AUTHOR, DATE
    INTEGER(KIND=PIOLONG) :: BEGININDX, ENDINDX
    
    
    myerr = PIOInfoObject(TOITYPE, DATATYPE, BEGININDX, ENDINDX, AUTHOR, DATE, file)
    if (myerr < 0_i8b) then
       print*,PIOErrMess(myerr)
       call fatal_error
    endif
    
    outstr = trim(toitype)
    
    return
  end function pio_objecttype
  !---------------------------
  function pio_datatype(file) result (outstr)

    character(len=MYLEN) :: outstr
    character(len=*), intent(in) :: file
    integer(piolong) :: myerr

    CHARACTER(LEN=MYLEN)  :: TOITYPE, DATATYPE, AUTHOR, DATE
    INTEGER(KIND=PIOLONG) :: BEGININDX, ENDINDX
    
    
    myerr = PIOInfoObject(TOITYPE, DATATYPE, BEGININDX, ENDINDX, AUTHOR, DATE, file)
    if (myerr < 0_i8b) then
       print*,PIOErrMess(myerr)
       call fatal_error
    endif
    
    outstr = trim(datatype)
    
    return
  end function pio_datatype
  
  !---------------------------
  subroutine check_pio_read(myerr, file, ndata, npix, code)

    integer(piolong), intent(in) :: myerr
    character(len=*), intent(in) :: file
    integer(i4b),     intent(inout) :: npix
    integer(i4b),     intent(in) :: ndata
    character(len=*), intent(in) :: code

    if (npix > 0) then
       !!     call assert_read_piofile(myerr, npix, file, code)
       if (myerr /= npix) then
          write(*,*) 'ERROR in '//code
          write(*,*) 'Inconsistent number of data for file'
          write(*,*) trim(file)
          write(*,*) 'Number of pixels detected :',myerr
          write(*,*) 'Npixtot should be :',npix
          CALL FATAL_ERROR
       endif
    endif

    if (myerr < 0_i8b) then
       print*,PIOErrMess(myerr)
       call fatal_error
    endif

    npix = min(myerr, ndata)

    return
  end subroutine check_pio_read

  ! ===========================================================================
  ! read arbitrary data type into single or double precision float array
  ! mypioreadmap(
  !            file, ! dmc object
  !            data, ! output vector, SP or DP
  !            npix, ! optional, number of elements expected to be read
  !            code) ! optional, code of calling routine
  ! ===========================================================================
  subroutine mypioreadmap_sp(file, data, npix, code)
  ! ===========================================================================
    integer, parameter :: KMAP = SP

    real(KMAP), intent(inout), dimension(1:) :: data
    character(len=*), intent(in) :: file
    integer(i4b),     intent(in), optional :: npix
    character(len=*), intent(in), optional :: code

    integer(piolong) :: myerr
    logical(1), pointer, dimension(:) :: pt_flag
    !integer(i1b), pointer, dimension(:) :: pt_flag
    integer(i1b), pointer, dimension(:) :: pt_byte ! never tested
    integer(i2b), pointer, dimension(:) :: pt_short ! never tested
    integer(i4b), pointer, dimension(:) :: pt_int ! never tested
    integer(i8b), pointer, dimension(:) :: pt_long ! never tested
    real(KMAP),   pointer, dimension(:) :: pt_float
    integer(i4b) :: my_ndata, my_npix
    character(len=MYLEN) :: my_code, datatype
    !------------------------------------------------------------------------

    my_code = ' '
    if (present(code)) my_code = code
    my_ndata = size(data)
    my_npix = 0_i4b
    if (present(npix)) my_npix = npix

    datatype = pio_datatype(file)
    select case (trim(datatype))
    case ('PIOFLAG')
       myerr = PIOreadmapobject(pt_flag,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_flag(1:my_npix)*1, kind=KMAP)
       myerr=PIOdeletemaptable(pt_flag)
    case ('PIOBYTE')
       myerr = PIOreadmapobject(pt_byte,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_byte(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_byte)
    case ('PIOSHORT')
       myerr = PIOreadmapobject(pt_short,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_short(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_short)
    case ('PIOINT')
       myerr = PIOreadmapobject(pt_int,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_int(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_int)
    case ('PIOLONG')
       myerr = PIOreadmapobject(pt_long,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_long(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_long)
    case ('PIOFLOAT','PIODOUBLE')
       myerr = PIOreadmapobject(pt_float,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_float(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_float)
    case default
       print*,'unexpected datatype: '//trim(datatype)
       call fatal_error
    end select

    return
  end subroutine mypioreadmap_sp

  ! ===========================================================================
  subroutine mypioreadmap_dp(file, data, npix, code)
  ! ===========================================================================
    integer, parameter :: KMAP = DP

    real(KMAP), intent(inout), dimension(1:) :: data
    character(len=*), intent(in) :: file
    integer(i4b),     intent(in), optional :: npix
    character(len=*), intent(in), optional :: code

    integer(piolong) :: myerr
    logical(1), pointer, dimension(:) :: pt_flag
    !integer(i1b), pointer, dimension(:) :: pt_flag
    integer(i1b), pointer, dimension(:) :: pt_byte ! never tested
    integer(i2b), pointer, dimension(:) :: pt_short ! never tested
    integer(i4b), pointer, dimension(:) :: pt_int ! never tested
    integer(i8b), pointer, dimension(:) :: pt_long ! never tested
    real(KMAP),   pointer, dimension(:) :: pt_float
    integer(i4b) :: my_ndata, my_npix
    character(len=MYLEN) :: my_code, datatype
    !------------------------------------------------------------------------

    my_code = ' '
    if (present(code)) my_code = code
    my_ndata = size(data)
    my_npix = 0_i4b
    if (present(npix)) my_npix = npix

    datatype = pio_datatype(file)
    select case (trim(datatype))
    case ('PIOFLAG')
       myerr = PIOreadmapobject(pt_flag,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_flag(1:my_npix)*1, kind=KMAP)
       myerr=PIOdeletemaptable(pt_flag)
    case ('PIOBYTE')
       myerr = PIOreadmapobject(pt_byte,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_byte(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_byte)
    case ('PIOSHORT')
       myerr = PIOreadmapobject(pt_short,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_short(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_short)
    case ('PIOINT')
       myerr = PIOreadmapobject(pt_int,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_int(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_int)
    case ('PIOLONG')
       myerr = PIOreadmapobject(pt_long,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_long(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_long)
    case ('PIOFLOAT','PIODOUBLE')
       myerr = PIOreadmapobject(pt_float,file,' ')
       call check_pio_read(myerr, file, my_ndata, my_npix, my_code)
       data(1:npix) = real(pt_float(1:my_npix), kind=KMAP)
       myerr=PIOdeletemaptable(pt_float)
    case default
       print*,'unexpected datatype: '//trim(datatype)
       call fatal_error
    end select

    return
  end subroutine mypioreadmap_dp

end module my_pio_routines

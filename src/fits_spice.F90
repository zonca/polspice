

module fits_spice

  use healpix_types
  use healpix_modules
  use spice_parameters, only: KMAP
  implicit none

  !integer(I4b), parameter :: KMAP = SP

  !-------------------------------------------------------------------------------
  ! parameters defined in Healpix fitstools, but kept private
  !-------------------------------------------------------------------------------
  real(kind=SP),     private, parameter :: s_bad_value = HPX_SBADVAL
  real(kind=DP),     private, parameter :: d_bad_value = HPX_DBADVAL
  integer(kind=I4B), private, parameter :: i_bad_value = -1637500000
  integer(I4B) ,     private, parameter :: nchunk_max  = 12000
  integer(I4B),      private, parameter :: MAXDIM_TOP  = 199 ! < 999

  interface f90ftgcv_
     module procedure f90ftgcve, f90ftgcvd
  end interface
  interface f90ftgky_
     module procedure f90ftgkye, f90ftgkyd
  end interface
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------


contains

  !-------------------------------------------------------------------------------
  ! routines defined in Healpix fitstools, but kept private
  !-------------------------------------------------------------------------------
  ! generic interface F90FTGCV_ for FITSIO's FTGCVE and FTGCVD
  !           reads data from BINTAB
  subroutine f90ftgcve(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(SP),     intent(out), dimension(0:) :: data
    real(SP),     intent(in)                 :: nullval
    call ftgcve(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcve
  subroutine f90ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(DP),     intent(out), dimension(0:) :: data
    real(DP),     intent(in)                 :: nullval
    call ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcvd
  !-------------------------------------------------------------------------------
  ! generic interface F90FTGKY_ for FITSIO's FTGKYE and FTGKYD
  !           reads a keyword
  subroutine f90ftgkye(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    real(SP),         intent(out) :: value
    call ftgkye(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkye
  subroutine f90ftgkyd(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    real(DP),         intent(out) :: value
    call ftgkyd(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkyd
  !-------------------------------------------------------------------------------
  !====================================================================
  subroutine get_clean_header(unit, header, filename, error, xalso, xonly)
    !====================================================================
    ! get_clean_header(unit, header, error [, xalso, xonly])
    ! gets the FITS header from unit, after excluding some default keywords 
    !  defined in def_excl
    ! if header in non void on input, its content will be concatenated with that
    !  of the FITS header
    ! if xalso is defined as a list of keywords, they are also excluded from the header
    ! if xonly is defined as a list of keywords, only those keywords are excluded from
    ! the header.
    ! xonly and xalso are exclusive
    !====================================================================
    INTEGER(I4B),                    intent(IN)           :: unit
    CHARACTER(LEN=*), DIMENSION(1:), INTENT(IN OUT)       :: header
    CHARACTER(LEN=*),                INTENT(IN)           :: filename
    INTEGER(I4B),                    intent(OUT)          :: error
    character(len=8), dimension(1:), intent(IN), optional :: xalso
    character(len=8), dimension(1:), intent(IN), optional :: xonly

    INTEGER(I4B) :: nlheader, status, i, n_excl
    CHARACTER(LEN=80) :: record
    CHARACTER(len=8), dimension(:), allocatable :: to_excl

    CHARACTER(len=8), dimension(1:21) :: def_excl
    !====================================================================

    ! keywords to be excluded by default from output header
    ! Note that TTYPE# and TUNIT# keywords are not excluded
    ! from the output header as they might be useful to the 
    ! calling routines
    def_excl=(/&
         & "SIMPLE  ","BITPIX  ","NAXIS   ",&
         & "NAXIS#  ","PCOUNT  ","GCOUNT  ",&
         & "EXTEND  ","ORIGIN  ","DATE*   ",&
         & "TFIELDS ","TFORM#  ",           & 
         & "TBCOL#  ","EXTNAME ","CTYPE#  ",&
         & "CRVAL#  ","CRPIX#  ","CDELT#  ",&
         & "XTENSION","INSTRUME","TELESCOP",&
         & "PDMTYPE "/)

    error = 0
    record=''

    if (present(xonly)) then 
       n_excl = size(xonly)
       allocate(to_excl(1:n_excl))
       to_excl = xonly

    else if (present(xalso)) then
       n_excl = size(xalso) + size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl(1:size(def_excl)) = def_excl
       to_excl(size(def_excl)+1:n_excl) = xalso

    else
       n_excl = size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl = def_excl
    endif

    nlheader=size(header)
    ! go to end of fortran header
    do i = 1, nlheader
       if (trim(header(i)) == "") exit
    enddo
    ! go to top of fits file header
    status=0
    call ftgrec(unit,0_i4b,record,status)
    ! read in all header lines except those excluded
    do
       call ftgnxk(unit,'*',1_i4b,to_excl,n_excl,record,status)
       if (status > 0) exit ! end of header
       if (i > nlheader) then
          write(unit=*,fmt="(a,i5,a)") &
               & " WARNING : The header in "//  &
               &    trim(filename)//" has more than ", &
               &  nlheader," lines."
          print*," It will be truncated."
          error = 1
          exit
       endif
       header(i)=record
       i=i+1
    enddo
    status=0

    return
  end subroutine get_clean_header
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  ! new routine(s)
  !-------------------------------------------------------------------------------
  !=======================================================================
  subroutine fix_polcconv(filename, map, polcconv_in, polcconv_out, header_in, header_out)
    !=======================================================================
    character(len=FILENAMELEN),       intent(in)            :: filename
    real(KMAP),     dimension(0:,1:), intent(inout)         :: map
    integer(I4B),                     intent(in)            :: polcconv_in
    integer(I4B),                     intent(out), optional :: polcconv_out
    character(len=80), dimension(1:), intent(in),  optional :: header_in
    character(len=80), dimension(1:), intent(out), optional :: header_out

    character(len=80) :: polcconv
    character(len=*), parameter :: primer_url = 'http://healpix.sf.net/pdf/intro.pdf'
    !------------------------------------------------------------------------

    select case (polcconv_in)
    case (0)
       print 9000,' The POLCCONV keyword was not found in '//trim(filename)
       print 9000,' COSMO (HEALPix default) will be assumed, and map is unchanged.'
       print 9000,' See HEALPix primer ('//primer_url//') for details.'
       
    case (1)
       continue
       ! print 9000,' POLCCONV=COSMO found in '//trim(filename)//'. Map is unchanged.'
          
    case (2)
       print 9000,' POLCCONV=IAU found in '//trim(filename)
       map(:,3) = - map(:,3)
       if (present(header_out)) then
          print 9000,' The sign of the U polarisation is changed in memory,'&
               & //' and the header is updated.'
          call add_card(header_out, 'POLCCONV', 'COSMO', &
               comment='Changed by input_map', update=.true.)
       else
          print 9000,' The sign of the U polarisation is changed in memory.'
       endif
       print 9000,' See HEALPix primer ('//primer_url//') for details.'
       
    case default
       if (present(header_in)) then
          call get_card(header_in,'POLCCONV',polcconv)
          print 9000,' POLCCONV='//trim(polcconv)//' found in '//trim(filename)
          print 9000,' It is neither COSMO nor IAU.   Aborting!'
       else
          print 9000,' The value of POLCCONV found in '//trim(filename)
          print 9000,' is neither COSMO nor IAU.   Aborting!'             
       endif
       print 9000,' See HEALPix primer ('//primer_url//') for details.'
       call fatal_error
    end select

    if (present(polcconv_out)) polcconv_out = 1

9000 format(a)
    return
  end subroutine fix_polcconv

  !=======================================================================
  subroutine read_fits_cutx(filename, pixel, cutmap, header, extno)
    !=======================================================================
    !
    ! routine to read FITS file with cut sky data : PIXEL, ?, ?, ?, ?, ...
    !
    ! read_fits_cutx(filename, pixel, cutmap, [header, extno])
    !=======================================================================

    character(len=*),                     intent(in)            :: filename
    integer(I4B),     dimension(0:),      intent(out)           :: pixel
    real(KMAP),       dimension(0:,1:),   intent(out)           :: cutmap
    character(len=*), dimension(1:),      intent(out), optional :: header
    integer(I4B),                         intent(in),  optional :: extno

    integer(I4B) :: obs_npix

    integer(I4B), parameter :: MAXDIM = MAXDIM_TOP !number of columns in the extension
    integer(I4B) :: blocksize, datacode
    integer(I4B) :: firstpix, frow, hdutype
    integer(I4B) :: naxis, nfound, nmove, npix, nrows, nmaps, nmr, i
    integer(I4B) :: readwrite, status, tfields, unit, varidat, width
    integer(I4B) :: repeat, repeat1, repeat2
    logical(LGT) :: anynull, extend
    character(len=20), dimension(MAXDIM) :: ttype, tform, tunit
    character(len=20) :: extname
    character(len=80) :: comment
    real(KMAP) :: nullval
    integer(i4b) :: npixtot
    !=======================================================================

    !npixtot = min(size_long(pixel),size_long(cutmap,1))
    npixtot = min(size(pixel),size(cutmap,1))
    nmaps   = size(cutmap,2)
    status=0

    unit = 150
    nfound = -1
    anynull = .false.

    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)
    if (naxis > 0) then ! there is an image
       print*,'an image was found in the FITS file '//filename
       print*,'... it is ignored.'
    endif

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) then 
       print*,'extension expected and not found in FITS file '//filename
       print*,'abort code'
       call fatal_error
    endif

    nmove = +1
    if (present(extno)) nmove = +1 + extno
    call ftmrhd(unit, nmove, hdutype, status)

    call assert (hdutype==2, 'this is not a binary table')

    ! reads all the FITS related keywords
    call ftghbn(unit, MAXDIM, &
         &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
         &        status)

    if (.not.  (trim(ttype(1)) == 'PIXEL' &
         & .or. trim(ttype(1)) == 'PIX'  ) ) call fatal_error('did not find PIXEL column in '//filename)
    

    if (present(header)) then 
       header = ""
       status = 0
       !call fitstools_mp_get_clean_header(unit, header, filename, status)
       call get_clean_header(unit, header, filename, status)
    endif

!     if (present(units)) then 
!        ! the second column contains the SIGNAL, for which we 
!        ! already read the units from the header
!        units = adjustl(tunit(2))
!     endif

    !        finds the bad data value
    call f90ftgky_(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (KMAP == SP) nullval = s_bad_value ! default value
       if (KMAP == DP) nullval = d_bad_value ! default value
       status = 0
    endif

    frow = 1
    firstpix = 1
    !        parse TFORM keyword to find out the length of the column vector
    repeat1 = 1
    repeat2 = 1
    call ftbnfm(tform(1), datacode, repeat1, width, status)
    if (tfields > 1) call ftbnfm(tform(2), datacode, repeat2, width, status)
    repeat = max(repeat1, repeat2)

    call ftgkyj(unit,'OBS_NPIX',obs_npix,comment,status)
    if (status == 202) then ! obs_npix not found
       obs_npix = nrows * repeat
       status = 0
    endif

    npix = min(npixtot, obs_npix)
    nmr  = min(nmaps,   tfields-1)
    call ftgcvj   (unit, 1_i4b, frow, firstpix, npix, i_bad_value, pixel(0), anynull, status)
    do i=1,nmr
       call f90ftgcv_(unit, i+1_i4b, frow, firstpix, npix, nullval, cutmap(0:npix-1,i), anynull, status)
    enddo

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_fits_cutx

end module fits_spice


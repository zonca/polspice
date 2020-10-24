!=======================================================================
subroutine convert_to_ascii(string,number)
!=======================================================================
  use healpix_types
  implicit none
  character(len=8) :: string,stringtmp
  integer(I4B) :: number

  if (number.gt.99999999.or.number.lt.0) then
     write(*,*) 'ERROR in convert_to_ascii'
     write(*,*) 'number > 99999999 or number < 0.'
     write(*,*) 'number =',number
     STOP
  endif
  write(string,'(I8.8)') number

end subroutine convert_to_ascii


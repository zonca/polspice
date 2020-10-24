program test_misc

! piece of code to test different things

  use misc_utils, only: fatal_error
  use test_misc_pie_io, only: init_parContent

  print*,"************************************************"
  print*,'Hello world'
  
!   print*,'call fatal_error'
!   call fatal_error
  
  print*,'call exit(-1)'
  call exit(-1)
  
  print*,"success"
  print*,"************************************************"
  stop

end program test_misc

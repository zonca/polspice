!ifort15hpxomp sets.F90 -c -o sets.o ; ifort15hpxomp -I. sets.o test_sets.F90 -o test_sets
!ifort15hpxomp -check all -traceback sets.F90 test_sets.F90 -o test_sets

program test_sets

  use healpix_types
  use sets, only: intersection, union
  implicit none
  integer(i4b), allocatable, dimension(:) :: list1, list2
  integer(i4b), dimension(:,:), pointer :: linter, lunion
  integer(i4b) :: i, j, n1, n2


  n1 = 10
  n2 = 5
  allocate(list1(0:n1-1))
  allocate(list2(0:n2-1))
  list1 = (/ (i*10, i=0,n1-1) /)
  list2 = (/ (i*10+10, i=0,n2-1) /)

!   n1 = 6
!   n2 = 5
!   allocate(list1(0:n1-1))
!   allocate(list2(0:n2-1))
!   list1 = (/ (i*10, i=n1-1,0,-1) /)
!   list2 = (/ (i*10, i=0,n2-1) /)

  print*,'List 1'
  print*,list1
  print*,'List 2'
  print*,list2

  linter => intersection(list1, list2)
  print*,'Intersection'
  print*,linter(:,0)
  print*,linter(:,1)
  print*,linter(:,2)

  lunion => union(list1, list2, sort=.true.)
  print*,'1 U 2'
  print*,lunion(:,0)
  print*,lunion(:,1)
  print*,lunion(:,2)

  lunion => union(list2, list1)
  print*,'2 U 1'
  print*,lunion(:,0)
  print*,lunion(:,1)
  print*,lunion(:,2)

end program test_sets

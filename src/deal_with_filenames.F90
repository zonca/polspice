!=======================================================================
subroutine deal_with_filenames
!=======================================================================
  use spice_common
  implicit none
  character(len=FILENAMELEN) :: filejunk_default,filejunk_input
  character(len=FILENAMELEN) :: ftmp_in, ftmp_def
  logical(lgt) :: Infile = .true., Outfile = .false., ltmp_def
  integer(i4b), dimension(1:3,1:2) :: listmaptype
  integer(i4b) :: jmap, jmmap, k
  integer(i4b) :: istokes
  
  filejunk_default=' ' ! unused

  listmaptype = 0 ! no map reading
  if (maxval(nlmaps(1:3,1)) == 0) then
     if (map_present)       listmaptype(1,1) = 1 ! mapfile only
     if (readQfile)         listmaptype(2,1) = 1
     if (readUfile)         listmaptype(3,1) = 1
     if (extramap1_present) listmaptype(1,1) = listmaptype(1,1)+2 ! mapfile + extramapfile
     if (readextraQfile1)   listmaptype(2,1) = listmaptype(2,1)+2
     if (readextraUfile1)   listmaptype(3,1) = listmaptype(3,1)+2
  else
     if (nlmaps(1,1) >0)    listmaptype(1,1) = 4 ! read from list, ignore actual (if any) mapfile_input and extramapfile
     if (nlmaps(2,1) >0)    listmaptype(2,1) = 4
     if (nlmaps(3,1) >0)    listmaptype(3,1) = 4
     map_present = .true. ! now indicates presence of any kind of map1
  endif
  if (maxval(nlmaps(1:3,2)) == 0) then
     if (map2_present)      listmaptype(1,2) = 1
     if (readQfile2)        listmaptype(2,2) = 1
     if (readUfile2)        listmaptype(3,2) = 1
     if (extramap2_present) listmaptype(1,2) = listmaptype(1,2)+2
     if (readextraQfile2)   listmaptype(2,2) = listmaptype(2,2)+2
     if (readextraUfile2)   listmaptype(3,2) = listmaptype(3,2)+2
  else
     if (nlmaps(1,2) >0)    listmaptype(1,2) = 4 ! read from list, ignore actual (if any) mapfile2_input and extramapfile2
     if (nlmaps(2,2) >0)    listmaptype(2,2) = 4
     if (nlmaps(3,2) >0)    listmaptype(3,2) = 4
     map2_present = .true. ! now indicates presence of any kind of map2
  endif

!   print*,map_present, extramap1_present,readextraQfile1,readextraUfile1
!   print*,map2_present,extramap2_present,readextraQfile2,readextraUfile2
!   print*,'listmaptype: '
!   print*,listmaptype(1:3,1)
!   print*,listmaptype(1:3,2)
  jmmap = 2

  ! First map and optional second or list of maps, for T(QU)
  do jmap=1,jmmap
     if (iand(listmaptype(1,jmap),1) == 1) then
        if (jmap == 1) then 
           ftmp_def = trim(mapfile_default)  ; ftmp_in = trim(mapfile_input)  ; ltmp_def = default_map_file
        endif
        if (jmap == 2) then 
           ftmp_def = trim(mapfile2_default) ; ftmp_in = trim(mapfile2_input) ; ltmp_def = default_map_file2
        endif
        call   get_generic_filename( &
             &    ftmp_def,          ftmp_in,                         listmapfile(1,1,jmap), &
             &         ltmp_def,        .true.,Infile,.false.,stringlength,mytype_junk)
        nlmaps(1,jmap) = 1
     endif
     if (iand(listmaptype(1,jmap),2)  == 2) then
        call    get_generic_filename(&
             &    filejunk_default,   extramapTQUfile_input(2,1,jmap), listmapfile(2,1,jmap), &
             &         .false.,         .true.,Infile,.false.,stringlength,mytype_junk)   
        nlmaps(1,jmap) = 2
     endif
     if (listmaptype(1,jmap)>3) then
        do k=1, nlmaps(1,jmap)
           call get_generic_filename(&
                & filejunk_default,   extramapTQUfile_input(k,1,jmap), listmapfile(k,1,jmap), &
                &      .false.,         .true.,Infile,.false.,stringlength,mytype_junk)  
        enddo
     endif
  enddo

#ifdef PIO
  do jmap=1,jmmap ! 1st (and optional second) map(s)
     do istokes=2,3 ! Q and U
        if (iand(listmaptype(istokes,jmap),1) == 1) then
           if (jmap==1 .and. istokes==2) ftmp_in = trim(mapQfile_input)
           if (jmap==1 .and. istokes==3) ftmp_in = trim(mapUfile_input)
           if (jmap==2 .and. istokes==2) ftmp_in = trim(mapQfile2_input)
           !if (jmap==2 .and. istokes==3) ftmp_in = trim(mapQfile2_input) ! bug corrected on 2015-10-12
           if (jmap==2 .and. istokes==3) ftmp_in = trim(mapUfile2_input)
           call  get_generic_filename(&
                &  filejunk_default, ftmp_in,                         listmapfile(1,istokes,jmap),&
                &          .false.,         .false.,Infile,.false.,stringlength,mytype_junk)
           nlmaps(istokes,jmap) = 1
        endif
        if (iand(listmaptype(istokes,jmap),2) == 2) then
           call  get_generic_filename(&
                &  filejunk_default, extramapTQUfile_input(2,istokes,jmap),  listmapfile(2,istokes,jmap), &
                &          .false.,         .false.,Infile,.false.,stringlength,mytype_junk)
           nlmaps(istokes,jmap) = 2
        endif
        if (listmaptype(istokes,jmap)>3) then
           do k=1, nlmaps(istokes,jmap)
              call get_generic_filename(&
                   &  filejunk_default, extramapTQUfile_input(k,istokes,jmap),listmapfile(k,istokes,jmap), &
                   &       .false.,         .false.,Infile,.false.,stringlength,mytype_junk)  
           enddo
        endif
     enddo
  enddo

#endif

  if (masks_present) then
     call get_generic_filename(maskfile_default, &
 &                             maskfile_input, &
 &                             maskfile,default_mask_file, &
 &                             .true.,Infile,.false.,stringlength,mytype_junk)

     call get_generic_filename(maskfile, & ! default:same as maskfile
 &                             maskfilep_input, & ! user provided
 &                             maskfilep,& ! final value
 &                             .not.maskfilep_toread, & ! take default if not deemed 'to read'
 &                             .false.,Infile,.false.,stringlength,mytype_junk)
     maskfilep_toread = (maskfile /= maskfilep) ! make sure filenames are indeed different
  endif

  if (masks2_present) then
     call get_generic_filename(maskfile2_default, &
 &                             maskfile2_input, &
 &                             maskfile2,default_mask_file2, &
 &                             .true.,Infile,.false.,stringlength,mytype_junk)

     call get_generic_filename(maskfile2, & ! default:same as maskfile2
 &                             maskfilep2_input, & ! user provided
 &                             maskfilep2,& ! final value
 &                             .not.maskfilep2_toread, & ! take default if not deemed 'to read'
 &                             .false.,Infile,.false.,stringlength,mytype_junk)
     maskfilep2_toread = (maskfile2 /= maskfilep2) ! make sure filenames are indeed different
  endif

  if (weights_present) then 
     call get_generic_filename(weightfile_default, &
 &                             weightfile_input, &
 &                             weightfile,default_weight_file, &
 &                             .true.,Infile,.false.,stringlength,mytype_junk)

     call get_generic_filename(weightfile, & ! default:same as weightfile
 &                             weightfilep_input, & ! user provided
 &                             weightfilep,& ! final value
 &                             .not.weightfilep_toread, & ! take default if not deemed 'to read'
 &                             .false.,Infile,.false.,stringlength,mytype_junk)
     weightfilep_toread = (weightfile /= weightfilep) ! make sure filenames are indeed different
  endif

  if (weights2_present) then 
     call get_generic_filename(weightfile2_default, &
 &                             weightfile2_input, &
 &                             weightfile2,default_weight_file2, &
 &                             .true.,Infile,.false.,stringlength,mytype_junk)

     call get_generic_filename(weightfile2, & ! default:same as weightfile2
 &                             weightfilep2_input, & ! user provided
 &                             weightfilep2,& ! final value
 &                             .not.weightfilep2_toread, & ! take default if not deemed 'to read'
 &                             .false.,Infile,.false.,stringlength,mytype_junk)
     weightfilep2_toread = (weightfile2 /= weightfilep2) ! make sure filenames are indeed different
  endif

  if (coroutput) call get_generic_filename(cor_file_name_default, &
 &                                         cor_file_name_input, &
 &                                         cor_file_name,default_cor_file, &
 &                                         .false.,Outfile,overwrite,stringlength,mytype_corr)

  if (cloutput) call get_generic_filename(cl_file_name_default, &
 &                                        cl_file_name_input, &
 &                                        cl_file_name,default_cl_file, &
 &                                        .false.,Outfile,overwrite,stringlength,mytype_cl)

  if (clmapoutput)  call get_generic_filename(clmap_file_name_default, &
 &                                            cl_outmap_file_input, &
 &                                            cl_outmap_file,default_clmapout_file, &
 &                                            .false.,Outfile,overwrite,stringlength,mytype_cl)

  if (clmapinput)  call get_generic_filename(clmap_file_name_default, &
 &                                           cl_inmap_file_input, &
 &                                           cl_inmap_file,default_clmapin_file, &
 &                                           .false.,Infile,.false.,stringlength,mytype_junk)


  if (clmaskoutput)  call get_generic_filename(clmask_file_name_default, &
 &                                            cl_outmask_file_input, &
 &                                            cl_outmask_file,default_clmaskout_file, &
 &                                            .false.,Outfile,overwrite,stringlength,mytype_cl)

  if (clmaskinput)  call get_generic_filename(clmask_file_name_default, &
 &                                            cl_inmask_file_input, &
 &                                            cl_inmask_file,default_clmaskin_file, &
 &                                            .false.,Infile,.false.,stringlength,mytype_junk)


  if (subtract_noise_cor) call get_generic_filename(noisecor_file_name_default, &
 &                                                   noisecor_file_name_input, &
 &                                                   noisecor_file_name, &
 &                                                   default_noisecor_file,.false., &
 &                                                   Infile,.false.,stringlength,mytype_corr)

  if (subtract_noise_cl) call get_generic_filename(noisecl_file_name_default, &
 &                                                  noisecl_file_name_input, &
 &                                                  noisecl_file_name, &
 &                                                  default_noisecl_file,.false., &
 &                                                  Infile,.false.,stringlength,mytype_cl)

  if (tenorminput) call get_generic_filename(tenorm_file_name_default, &
 &                                           tenorm_in_file_input, &
 &                                           tenorm_in_file, &
 &                                           default_tenormin_file,.true., &
 &                                           Infile,.false.,stringlength,mytype_cl)

  if (tenormoutput) call get_generic_filename(tenorm_file_name_default, &
 &                                            tenorm_out_file_input, &
 &                                            tenorm_out_file, &
 &                                            default_tenormout_file,.true., &
 &                                            Outfile,overwrite,stringlength,mytype_cl)

  if (windowinput) call get_generic_filename(window_file_name_default, &
 &                                           window_in_file_input, &
 &                                           window_in_file, &
 &                                           default_windowin_file,.true., &
 &                                           Infile,.false.,stringlength,mytype_table)

  if (windowoutput) call get_generic_filename(window_file_name_default, &
 &                                            window_out_file_input, &
 &                                            window_out_file, &
 &                                            default_windowout_file,.true., &
 &                                            Outfile,overwrite,stringlength,mytype_table)

  if (kernelsoutput) call get_generic_filename(kernels_file_name_default, &
 &                                            kernels_out_file_input, &
 &                                            kernels_out_file, &
 &                                            default_kernelsout_file,.true., &
 &                                            Outfile,overwrite,stringlength,mytype_table3)

  if (do_cov) call get_generic_filename(cov_file_name_default, &
 &                                            cov_out_file_input, &
 &                                            cov_out_file, &
 &                                            default_covout_file,.true., &
 &                                            Outfile,overwrite,stringlength,mytype_table3)

  if (beam_present) then
     filejunk_input=trim(beam_file)
     call get_generic_filename(filejunk_default, &
 &                             filejunk_input, &
 &                             beam_file,.false.,.false.,Infile,.false.,stringlength,mytype_junk)                              
  endif

  if (beam2_present) then
     filejunk_input=trim(beam_file2)
     call get_generic_filename(filejunk_default, &
 &                             filejunk_input, &
 &                             beam_file2,.false.,.false.,Infile,.false.,stringlength,mytype_junk)                              
  endif
end subroutine deal_with_filenames

!=======================================================================
subroutine get_generic_filename(file_name_default, &
 &                              file_name_input, &
 &                              file_name, default, &
 &                              fits,in,overwrite,stringlength,mytype)
!=======================================================================
! checks that file name is valid and that file
! - exists if trying to read it (in == .true.)
! - can be written/overwritten if trying to do so (overwrite==.true.)
! if default=.true., force file_name to file_name_default (adding a '.fits' suffix is fits==.true.) otherwise use file_name_input
!=======================================================================
!
  use healpix_types
  use misc_utils, only: fatal_error, file_present
#ifdef PIO
  use piolib
  use spice_common, only : out_cl_grp ! EH
  use spice_parameters, only : mytype_vect, mytype_cl, mytype_corr, mytype_table, &
       & mytype_vect_tuple, mytype_cl_tuple, mytype_corr_tuple, &
       & know_out_cl_grp, mytype_table3
#endif
  implicit none
  
  character(len=*), intent(IN) :: file_name_default
  character(len=*), intent(IN) :: file_name_input
  character(len=FILENAMELEN), intent(OUT) :: file_name
  logical(LGT),     intent(IN) :: default,fits,in,overwrite
  integer(I4B),     intent(IN) :: stringlength
  integer(I4B),     intent(IN) :: mytype

  integer(I4B) :: leng,lengsub
  logical :: exist, exist_extended, create_and_delete = .false.

#ifdef PIO
  integer(I8B) :: myerr,myerr2,mygroup
  character(len=DMCPIOSTRINGMAXLEN) :: grpname
!  character*128 :: grpname
#endif

  !!!!print*,'in get_generic_filename',trim(file_name_default),trim(file_name_input),default,fits,in,overwrite,mytype
  if (default) then
     leng=len(trim(file_name_default))
     if (fits) then
        lengsub=5
     else
        lengsub=0
     endif
     if (leng > stringlength-lengsub) then
        write(*,*) 'ERROR in get_generic_filename'
        write(*,*) 'file name too long :'
        write(*,*) trim(file_name_default)
        CALL FATAL_ERROR
     endif
     file_name=trim(file_name_default)
#ifndef PIO
     if (fits)  file_name=trim(file_name_default)//'.fits'
#endif
  else
     leng=len(trim(file_name_input))
     if (leng > stringlength) then
        write(*,*) 'ERROR in get_generic_filename'
        write(*,*) 'file name too long :'
        write(*,*) trim(file_name_input)
        CALL FATAL_ERROR
     endif
     file_name=trim(file_name_input)
  endif
     
#ifdef PIO
  if (in) then ! inputting file
     myerr=PIOCheckObject(file_name)
     if (myerr /= 0) then
        write(*,*) 'ERROR ',myerr,' in get_generic_filename for input file'
        write(*,*) trim(file_name)
        write(*,*) PIOerrmess(myerr)
        CALL FATAL_ERROR
     endif
  else ! outputting file
     myerr=PIOCheckObject(file_name)
     if (myerr == 0 .and. .not. overwrite) then
        write(*,*) 'ERROR in get_generic_filename'
        write(*,*) 'This object already exists :'
        write(*,*) trim(file_name)
        CALL FATAL_ERROR
     elseif (myerr /= 0) then ! file does not exist yet
        ! particular case of temporary vector file
        ! issue warning, do not attempt to create nor open group
        if ((mytype==mytype_cl .or. mytype==mytype_vect .or. mytype==mytype_corr) &
             &  .and. PIOCheckTemporary(file_name) == 0 &
             &  .and. know_out_cl_grp) then
           grpname = trim(out_cl_grp)
           write(*,*) 'WARNING: Temporary file, user provided group '//trim(grpname)//' will be used' 
           return
        endif
        myerr=PIOgetgrpname(grpname,file_name)
        select case (mytype)
        case (mytype_vect, mytype_vect_tuple) ! generic: vect
           myerr   = PIOCreateVectGrp(grpname)
           mygroup = PIOOpenVectGrp(grpname,'w')
           myerr2  = PIOCloseVectGrp(mygroup)

           if (create_and_delete) myerr=PIOcreatevectobject(file_name,'PIODOUBLE')
           if (create_and_delete) myerr2=PIOdeleteobject(file_name)

        case (mytype_cl, mytype_cl_tuple, mytype_corr, mytype_corr_tuple) ! particular: cl and corr
           myerr   = PIOCreateClGrp(grpname)
           mygroup = PIOOpenClGrp(grpname,'w')
           myerr2  = PIOCloseClGrp(mygroup)

           if (create_and_delete) myerr=PIOcreateClobject(file_name,'PIODOUBLE')
           if (create_and_delete) myerr2=PIOdeleteobject(file_name)

        case (mytype_table)
           myerr   = PIOCreateTab2dGrp(grpname,1_I8B)
           mygroup = PIOOpenTab2dGrp(grpname,'w')
           myerr2  = PIOCloseTab2dGrp(mygroup)

           if (create_and_delete) myerr=PIOcreatetab2dobject(file_name,'PIODOUBLE')
           if (create_and_delete) myerr2=PIOdeleteobject(file_name)

        case (mytype_table3)
!            myerr   = PIOCreateTab3dGrp(grpname,1_I8B,1_I8B)
!            mygroup = PIOOpenTab3dGrp(grpname,'w')
!            myerr2  = PIOCloseTab3dGrp(mygroup)

!            if (create_and_delete) myerr=PIOcreatetab3dobject(file_name,'PIODOUBLE')
!            if (create_and_delete) myerr2=PIOdeleteobject(file_name)
           myerr2 = 0
           mygroup = 0
        end select
        if (mygroup < 0 .or. myerr2 /= 0) then
           write(*,*) 'ERROR in get_generic_filename for output file '//trim(file_name)
           write(*,*) 'Group name :  '//trim(grpname)
           write(*,*) 'create group: '//PIOerrmess(myerr)
           write(*,*) 'open group:   '//PIOerrmess(mygroup)
           write(*,*) 'close group:  '//PIOerrmess(myerr2)
           CALL FATAL_ERROR
        endif
     endif
  endif
#else
  inquire(file=file_name,exist=exist)
  exist_extended = file_present(file_name)
  if (in .and. (.not.exist_extended)) then
     write(*,*) 'ERROR in get_generic_filename'
     write(*,*) 'INPUT file '//trim(file_name)
     write(*,*) 'does not exist.'
     print*,'default: '//trim(file_name_default), &
          & ', input: '//trim(file_name_input)
     CALL FATAL_ERROR
  elseif ((.not.in) .and. exist .and. (.not. overwrite)) then
     write(*,*) 'ERROR in get_generic_filename'
     write(*,*) 'OUTPUT file '//trim(file_name)
     write(*,*) 'already exists.'
     CALL FATAL_ERROR
  endif
#endif

end subroutine get_generic_filename


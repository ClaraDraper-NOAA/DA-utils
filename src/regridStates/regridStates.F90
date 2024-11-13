 program regridStates
! Program to re-grid a list of FV3 variables.
! Intended for use in DA applications (regridding of restarts for recentering, regridding increments).
!
! Clara Draper, and George Gayno  Aug, 2024. 

! TO DO - add option to pre-compute, and read in regridding route.

 use mpi_f08
 use esmf

 use grids_IO, only     : setup_grid, &
                          write_from_fields, &
                          read_into_fields, &
                          n_tiles, &
                          grid_setup_type

 use utilities, only    : error_handler

 implicit none

 integer, parameter             :: max_vars = 10  ! increase if wish to specify more variables
 
 ! namelist inputs
 character(len=100)             :: dir_fix
 integer                        :: ires_in, jres_in, ires_out, jres_out
 character(len=100)             :: dir_in, dir_out, dir_out_mask
 character(len=100)             :: fname_in, fname_out, fname_out_mask
 character(len=15)              :: variable_list(max_vars)
 character(len=10)              :: mask_type
 character(len=7)               :: gridtype_in, gridtype_out
 integer                        :: n_vars
 real(esmf_kind_r8)             :: missing_value ! value given to unmapped cells in the output grid

 type(grid_setup_type)          :: grid_setup_in, grid_setup_out

 integer                        :: ierr, localpet, npets
 integer                        :: v, SRCTERM
 
 type(esmf_vm)                  :: vm
 type(esmf_grid)                :: grid_in, grid_out
 type(esmf_field), allocatable  :: fields_in(:)
 type(esmf_field), allocatable  :: fields_out(:) 
 type(esmf_routehandle)         :: regrid_route
 real(esmf_kind_r8), pointer    :: ptr_out(:,:)



 real :: t1, t2, t3, t4

 ! see README for details of namelist variables.
 namelist /config/ fname_in, dir_in, gridtype_in, &
                   fname_out, dir_out, gridtype_out, &
                   fname_out_mask, dir_out_mask, &
                   n_vars, variable_list, missing_value, mask_type, &
                   ires_in, jres_in, ires_out, jres_out, &
                   dir_fix ! only needed for fv3* grids

!-------------------------------------------------------------------------
! INITIALIZE
!-------------------------------------------------------------------------

 call cpu_time(t1)

 call mpi_init(ierr)

 call ESMF_Initialize(rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("INITIALIZING ESMF", ierr)

 call ESMF_VMGetGlobal(vm, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGetGlobal", ierr)

 call ESMF_VMGet(vm, localPet=localpet, petCount=npets, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN VMGet", ierr)

!-------------------------------------------------------------------------
! RUN
!-------------------------------------------------------------------------

 print*,'** pets: local, total: ',localpet, npets

 ! checks 

 if (mod(npets,n_tiles) /= 0) then
   call error_handler("must run with a task count that is a multiple of 6", 1)
 endif

!------------------------
! read in namelist

 missing_value=-999. ! set defualt

 open(41, file='regrid.nml', iostat=ierr)
 if (ierr /= 0) call error_handler("OPENING regrid NAMELIST.", ierr)
 read(41, nml=config, iostat=ierr)
 if (ierr /= 0) call error_handler("READING regrid NAMELIST.", ierr)
 close (41)

  
 grid_setup_in%descriptor = gridtype_in

 ! if in
     grid_setup_in%dir_mask = dir_in ! for fv3 and gauss use input file for mask
     grid_setup_in%fname_mask = fname_in 
 grid_setup_in%dir = dir_in ! for fv3 and gauss use input file for mask
 grid_setup_in%fname = fname_in 
 ! if gau_inc
 !grid_setup_in%dir_coord = dir_in
 grid_setup_in%dir_coord = dir_fix
 !grid_setup_in%fname_coord = fname_in

 grid_setup_out%descriptor = gridtype_out
 ! if out
     grid_setup_out%dir_mask = dir_out_mask
     grid_setup_out%fname_mask = fname_out_mask
 grid_setup_out%dir = dir_out
 grid_setup_out%fname = fname_out
 ! if fv3 (will construct mosaic fname later)
 grid_setup_out%dir_coord = dir_fix

!------------------------
! Create esmf grid objects for input and output grids, and add land masks

! TO DO - can we make the number of tasks more flexible for fv3

 if (localpet==0) print*,'** Setting up grids'
 call setup_grid(localpet, npets, grid_setup_in, & 
                  trim(mask_type), grid_in, &
                  ires_in, jres_in)
 
 call setup_grid(localpet, npets, grid_setup_out, & 
                  trim(mask_type), grid_out, &
                  ires_out, jres_out)

!------------------------
! Create input and output fields

! TODO - think about 3-D fields

 if (localpet==0) print*,'** Creating/Reading fields'

! input
 allocate(fields_in(n_vars)) 

 do v = 1, n_vars

    fields_in(v)  = ESMF_FieldCreate(grid_in, &
                        typekind=ESMF_TYPEKIND_R8, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, &
                        name="input for regridding", &
                        rc=ierr)

    if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
       call error_handler("in FieldCreate "//trim(variable_list(v)), ierr)
 end do

! output
 allocate(fields_out(n_vars)) 

 do v = 1, n_vars

     fields_out(v)  = ESMF_FieldCreate(grid_out, &
                                       typekind=ESMF_TYPEKIND_R8, &
                                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                                       name="output of regridding", &
                                       rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("in FieldCreate, field_out", ierr)


     ! set the default output value (for non-mapped cells) 
     call ESMF_FieldGet(fields_out(v), &
                        farrayPtr=ptr_out, &
                        rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
     call error_handler("IN FieldGet", ierr)

     ptr_out=missing_value

 enddo

!------------------------
! read data into input fields

 call read_into_fields(localpet, ires_in, jres_in, trim(grid_setup_in%fname), &
                         trim(grid_setup_in%dir), grid_setup_in, n_vars,      &
                         variable_list(1:n_vars), fields_in) 

 call cpu_time(t2)
!------------------------
! regrid the input fields to the output grid

 if (localpet==0) print*,'** Performing regridding'

 SRCTERM=1
 ! get regriding route for a field (only uses the grid info in the field)
 ! to turn off masking, remove [src/dstMaskVales] argumemnts
 call ESMF_FieldRegridStore(srcField=fields_in(1), srcMaskValues=(/0/), &
                            dstField=fields_out(1), dstMaskValues=(/0/), &
                            ! allow unmapped grid cells, without returning error
                            unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
                            polemethod=ESMF_POLEMETHOD_ALLAVG, &
                            ! fill un-mapped grid cells with a neighbour
                            extrapMethod=ESMF_EXTRAPMETHOD_CREEP, & 
                            ! number of "levels" of neighbours to search for a value
                            extrapNumLevels=2, &
                            ! needed for reproducibility
                            ! (combined with ESMF_TERMORDER_SRCSEQ below)
                            srctermprocessing=SRCTERM, & 
                            routehandle=regrid_route, &
                            ! use bilinear interp (slightly better results than PATCH)
                            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridStore", ierr)

! do the re-gridding

 call cpu_time(t3)

 do v=1, n_vars
     call ESMF_FieldRegrid(fields_in(v), &
                           fields_out(v), &
                           routehandle=regrid_route, &
                           zeroregion=ESMF_REGION_SELECT, & ! initialize output with missing_value
                           termorderflag=ESMF_TERMORDER_SRCSEQ, rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
      call error_handler("IN FieldRegrid", ierr)
 enddo

! TO-DO: terrain-correct temperatures (all layers?)

! write out fields on destination grid

 if (localpet==0) print*,'** Writing out regridded fields'

 call write_from_fields(localpet, ires_out, jres_out, trim(grid_setup_out%fname),    &
                          trim(grid_setup_out%dir), n_vars, variable_list(1:n_vars), &
                          fields_out)


! clean up 

 call ESMF_FieldRegridRelease(routehandle=regrid_route, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN FieldRegridRelease", ierr)
 
 do v = 1, n_vars
     call ESMF_FieldDestroy(fields_in(v),rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
            call error_handler("DESTROYING FIELD", ierr)

     call ESMF_FieldDestroy(fields_out(v),rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("DESTROYING FIELD", ierr)
 enddo

 call ESMF_GridDestroy(grid_in,rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("DESTROYING GRID", ierr)

 call ESMF_GridDestroy(grid_out,rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("DESTROYING GRID", ierr)
 
!-------------------------------------------------------------------------
! FINALIZE
!-------------------------------------------------------------------------

 call ESMF_finalize(endflag=ESMF_END_KEEPMPI, rc=ierr)
 call mpi_finalize(ierr)

 call cpu_time(t4)
 if (localpet==0) print*, '** time in tile2tile', t4 - t1
 if (localpet==0) print*, '** time in RegridStore', t3 - t2

 print*,"** DONE.", localpet

 end program regridStates

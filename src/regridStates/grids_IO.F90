 module grids_IO
 ! ESMF grid-specific routines for GFS regridding program, including IO.
 ! 
 ! Clara Draper, Aug 2024.

 use esmf
 use netcdf
 use ESMF_LogPublicMod
 use utilities, only     : error_handler, netcdf_err

 implicit none

 private

 integer, public, parameter  :: n_tiles=6 ! number tiles in fv3 grid
 integer, public, parameter  :: vtype_water=0, & ! TO DO - which veg classification is this?
                                vtype_landice=15 ! used for soil mask
 type, public  :: grid_setup_type
        character(7)   :: descriptor
        character(100) :: fname
        character(100) :: dir
        character(100) :: fname_mask
        character(100) :: dir_mask
        character(100) :: fname_coord
        character(100) :: dir_coord
        integer        :: ires
        integer        :: jres
 end type
 
 public :: setup_grid, &
           read_into_fields, &
           write_from_fields

 contains

!-----------------------------------
! Create ESMF grid objects, with mask if requested

 subroutine setup_grid(localpet, npets, grid_setup,  & 
                       mask_type, mod_grid )

 implicit none

 ! INTENT IN
 type(grid_setup_type), intent(in)    :: grid_setup
 character(*), intent(in)       :: mask_type
 integer, intent(in)            :: localpet, npets

 ! INTENT OUT 
 type(esmf_grid)                :: mod_grid

 ! LOCAL
 type(esmf_field)               :: mask_field(1)
 real(esmf_kind_r8), pointer    :: ptr_maskvar(:,:)
 integer(esmf_kind_i4), pointer :: ptr_mask(:,:)

 integer                        :: ierr, ncid, tile
 character(len=15)              :: mask_variable(1)

!--------------------------
! Create grid object, and set up pet distribution

 select case (grid_setup%descriptor)
 case ('fv3_rst')
     call create_grid_fv3(grid_setup%ires, trim(grid_setup%dir_coord), npets, localpet ,mod_grid)
     mask_variable(1) = 'vtype          '
 case ('gau_inc')
     call create_grid_gauss(grid_setup, npets, localpet,  mod_grid)
     mask_variable(1) = 'soilsnow_mask  '
 case default 
     call error_handler("unknown grid_setup%descriptor in setup_grid", 1)
 end select

!--------------------------
! Calculate and add the mask

 if (mask_type=="soil") then 
 ! mask out ocean and glaciers, using vegetation class

     mask_field(1) = ESMF_FieldCreate(mod_grid, &
                                       typekind=ESMF_TYPEKIND_R8, &
                                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                                       name="input variable for mask", &
                                       rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate, mask_variable", ierr)

     call read_into_fields(localpet, grid_setup%ires, grid_setup%jres, trim(grid_setup%fname_mask), &
                             trim(grid_setup%dir_mask), grid_setup, 1, &
                             mask_variable(1), mask_field(1))

    ! get pointer to vegtype
     call ESMF_FieldGet(mask_field(1), &
                        farrayPtr=ptr_maskvar, &
                        rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGet", ierr)

    ! create and get pointer to the mask
     call ESMF_GridAddItem(mod_grid, &
                           itemflag=ESMF_GRIDITEM_MASK, &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("in GridAddItem mask", ierr)

     call ESMF_GridGetItem(mod_grid, & 
                           itemflag=ESMF_GRIDITEM_MASK, &
                           farrayPtr=ptr_mask, &
                           rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("in GridGetItem mask", ierr)

    ! calculate the mask
     ptr_mask = 1 ! initialize land everywhere
     if (mask_variable(1) == 'vtype          ') then
         where (nint(ptr_maskvar) == vtype_water )   ptr_mask = 0 ! exclude water
         where (nint(ptr_maskvar) == vtype_landice ) ptr_mask = 0 ! exclude glaciers
     end if

    ! destroy veg type field
     call ESMF_FieldDestroy(mask_field(1),rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("DESTROYING FIELD", ierr)

  end if ! mask = soil 

 end subroutine setup_grid

 ! read variables from fv3 netcdf restart file into ESMF Fields
 subroutine read_into_fields(localpet, i_dim, j_dim , fname_read, dir_read, &
                               grid_setup, n_vars, variable_list, fields) 
 
 implicit none 

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_read
 character(*), intent(in)        :: dir_read
 type(grid_setup_type), intent(in)        :: grid_setup

 character(len=15), dimension(n_vars), intent(in)   :: variable_list
 
 ! INTENT OUT 
 type(esmf_field), dimension(n_vars), intent(inout) :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, v
 integer                         :: n_files
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_in(:,:,:)

 if (localpet==0) then
     allocate(array_in(n_vars,i_dim, j_dim))
     allocate(array2D(i_dim, j_dim))
 else 
     allocate(array_in(0,0,0))
     allocate(array2D(0,0))
 end if

 select case (grid_setup%descriptor)
 case ('fv3_rst')
     n_files=n_tiles
 case ('gau_inc') 
     n_files=1
 case default
     call error_handler("unknown grid_setup%descriptor in read into fields", 1)
 end select
     
 do tt = 1, n_files

      ! read from restart
      if (localpet == 0) then

         if ( n_files > 1) then
             write(tchar,'(i1)') tt
             fname = dir_read//"/"//fname_read//".tile"//tchar//".nc"
         else
             fname = dir_read//"/"//fname_read
         endif

         !print *, 'Reading ', trim(fname)

         ierr=nf90_open(trim(fname),NF90_NOWRITE,ncid)
         call netcdf_err(ierr, 'opening: '//trim(fname) )

         do v =1, n_vars
             !print *, 'Reading ', trim(variable_list(v))
             ierr=nf90_inq_varid(ncid, trim(variable_list(v)), id_var)
             call netcdf_err(ierr, 'reading variable id' )

             ierr=nf90_get_var(ncid, id_var, array_in(v,:,:))
             call netcdf_err(ierr, 'reading variable' )
         enddo
         ierr = nf90_close(ncid)

      endif

      ! scatter
      do v =1, n_vars
          array2D=array_in(v,:,:) ! scatter misbehaves if given indexed 3D array.
          call ESMF_FieldScatter(fields(v), array2D, rootpet=0, tile=tt, rc=ierr)
          if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldScatter", ierr)
      enddo

 enddo
 
 ! clean up
 deallocate(array_in)

 end subroutine read_into_fields

 ! read lat and lon from history file (Guassian grid), and put into 2D fields
 ! needed for setting up the Gauusian grid.

 subroutine lonlat_read_into_fields(localpet, i_dim, j_dim, fname_read, dir_read, &
                               grid_setup, gauss_grid, n_vars, lonvar_list, &
                               latvar_list, lon_fields, lat_fields)
 
 implicit none 

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_read
 character(*), intent(in)        :: dir_read
 type(grid_setup_type), intent(in)        :: grid_setup
 type(esmf_grid), intent(in)     :: gauss_grid

 character(len=15), dimension(n_vars), intent(in)   :: lonvar_list
 character(len=15), dimension(n_vars), intent(in)   :: latvar_list
 
 ! INTENT OUT 
 type(esmf_field), dimension(n_vars), intent(out) :: lon_fields
 type(esmf_field), dimension(n_vars), intent(out) :: lat_fields

 ! LOCAL
 integer                         :: id_var, ncid, ierr
 integer                         :: i,j,v
 character(len=500)              :: fname
 real(esmf_kind_r8)              :: vec_lat(i_dim*j_dim), vec_lon(i_dim*j_dim)
 real(esmf_kind_r8), allocatable :: array_lat(:,:,:), array_lon(:,:,:), array2d(:,:)
 character(len=30)               :: vname

 ! Create the fields 
 do v=1,n_vars
        vname="gauss_grid_" // trim(lonvar_list(v)) 
        lon_fields(v) = ESMF_FieldCreate(gauss_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name=vname, rc=ierr)
        if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate", ierr)
 
        vname="gauss_grid_" // trim(latvar_list(v))
        lat_fields(v) = ESMF_FieldCreate(gauss_grid, &
                                   typekind=ESMF_TYPEKIND_R8, &
                                   staggerloc=ESMF_STAGGERLOC_CENTER, &
                                   name=vname, rc=ierr)
        if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate", ierr)
 enddo
 
 if (localpet==0) then
     allocate(array_lat(n_vars, i_dim, j_dim))
     allocate(array_lon(n_vars, i_dim, j_dim))
     allocate(array2d( i_dim, j_dim))
 else 
     allocate(array_lat(n_vars,0,0))
     allocate(array_lon(n_vars,0,0))
     allocate(array2d(0,0))
 endif

 ! read from restart
 if (localpet == 0) then

     fname = dir_read//"/"//fname_read

     print *, 'Reading ', trim(fname)

     ierr=nf90_open(trim(fname),NF90_NOWRITE,ncid)
     call netcdf_err(ierr, 'opening: '//trim(fname) )

     do v=1, n_vars

         print *, 'Reading ', trim(lonvar_list(v))
         ierr=nf90_inq_varid(ncid, trim(lonvar_list(v)), id_var)
         call netcdf_err(ierr, 'reading variable id' )

         ierr=nf90_get_var(ncid, id_var, vec_lon)
         call netcdf_err(ierr, 'reading variable' )

         array_lon(v,:,:) = reshape(vec_lon,(/i_dim,j_dim/))

         !do j = 1, j_dim
         !      array_lon(v,:,j) =vec_lon
         !enddo
         print *, 'CSD lons', array_lon(v,1,1), array_lon(v,i_dim,1)
         print *, 'CSD lons', array_lon(v,1,j_dim), array_lon(v,i_dim,j_dim)

         print *, 'Reading ', trim(latvar_list(v))
         ierr=nf90_inq_varid(ncid, trim(latvar_list(v)), id_var)
         call netcdf_err(ierr, 'reading variable id' )

         ierr=nf90_get_var(ncid, id_var, vec_lat)
         call netcdf_err(ierr, 'reading variable' )
         ierr = nf90_close(ncid)

         array_lat(v,:,:) = reshape(vec_lat,(/i_dim,j_dim/))
         !do i = 1, i_dim
         !      array_lat(v,i,:) =vec_lat
         !enddo
         print *, 'CSD lats', array_lat(v,1,1), array_lat(v,i_dim,1)
         print *, 'CSD lats', array_lat(v,1,j_dim), array_lat(v,i_dim,j_dim)
     enddo
     ierr = nf90_close(ncid)
  endif

  do v =1,n_vars

      ! scatter longitudes
      if (localpet==0) print *, 'scattering lon', lonvar_list(v)
      array2d=array_lon(v,:,:)
      call ESMF_FieldScatter(lon_fields(v), array2d, rootpet=0,  rc=ierr)
      if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldScatter", ierr)
     
      ! scatter latitude
      if (localpet==0) print *, 'scattering lat', latvar_list(v)
      array2d=array_lat(v,:,:)
      call ESMF_FieldScatter(lat_fields(v), array2d, rootpet=0,  rc=ierr)
      if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
         call error_handler("IN FieldScatter", ierr)
 enddo

 ! clean up
 deallocate(array_lat)
 deallocate(array_lon)
 deallocate(array2d)

 end subroutine lonlat_read_into_fields
! write variables from ESMF Fields into netcdf restart-like file 

 subroutine write_from_fields(localpet, i_dim, j_dim , fname_out, dir_out, &
                                n_vars, variable_list, fields) 

 implicit none 

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_out
 character(*), intent(in)        :: dir_out
 character(15), dimension(n_vars), intent(in)     :: variable_list
 type(esmf_field), dimension(n_vars), intent(in)  :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, &
                                   id_x, id_y, v
 character(len=1)                :: tchar
 character(len=500)              :: fname
 real(esmf_kind_r8), allocatable :: array2D(:,:)
 real(esmf_kind_r8), allocatable :: array_out(:,:,:)

 do v = 1, n_vars
        if (localpet == 0)  print *, 'Writing ', trim(variable_list(v)), ' into field'
 enddo

 if (localpet==0) then
     allocate(array_out(n_vars, i_dim, j_dim))
     allocate(array2D(i_dim, j_dim))
 else 
     allocate(array_out(0,0,0))
     allocate(array2D(0,0))
 end if

 do tt = 1, n_tiles

      ! fetch data (all PETs)
      do  v=1, n_vars
          call ESMF_FieldGather(fields(v), array2D, rootPet=0, tile=tt, rc=ierr)
          if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
             call error_handler("IN FieldGather", ierr)
          array_out(v,:,:) = array2D
      enddo

      ! write to netcdf 
      if (localpet == 0) then

         ! open file, set dimensions
         write(tchar,'(i1)') tt
         fname = dir_out//"/"//fname_out//".tile"//tchar//".nc"

         ierr = nf90_create(trim(fname), NF90_NETCDF4, ncid)
         call netcdf_err(ierr, 'creating file='//trim(fname) )

         ierr = nf90_def_dim(ncid, 'xaxis_1', i_dim, id_x)
         call netcdf_err(ierr, 'defining xaxis dimension' )

         ierr = nf90_def_dim(ncid, 'yaxis_1', j_dim, id_y)
         call netcdf_err(ierr, 'defining yaxis dimension' )

         do v=1, n_vars

             ierr = nf90_def_var(ncid, trim(variable_list(v)), NF90_DOUBLE, (/id_x, id_y/), id_var)
             call netcdf_err(ierr, 'defining '//variable_list(v) )
              
             ierr = nf90_put_var( ncid, id_var, array_out(v,:,:) ) 
             call netcdf_err(ierr, 'writing '//variable_list(v) ) 

         enddo

         ierr = nf90_close(ncid)

      endif

 enddo
 
 ! clean up
 deallocate(array_out)

 end subroutine write_from_fields

 ! subroutine to create grid object for fv3 grid
 ! also sets distribution across procs 

 subroutine create_grid_fv3(res_atm, dir_fix, npets, localpet, fv3_grid)

! INTENT IN
 integer, intent(in)    :: npets, localpet
 integer, intent(in)    :: res_atm
 character(*), intent(in)       :: dir_fix

 ! INTENT OUT 
 type(esmf_grid)                :: fv3_grid

 integer                :: ierr, extra, tile 
 integer                :: decomptile(2,n_tiles)
 
 character(len=5)       :: rchar
 character(len=200)     :: fname, dir_fix_res

 if (localpet == 0) print*," creating fv3 grid for ", res_atm 

! pet distribution
 extra = npets / n_tiles
 do tile = 1, n_tiles
   decomptile(:,tile)=(/1,extra/)
 enddo

 ! mosaic file 
 write(rchar,'(i3)') res_atm
 dir_fix_res = dir_fix//"/C"//trim(adjustl(rchar))//"/"
 fname = trim(dir_fix_res)//"/C"//trim(adjustl(rchar))// "_mosaic.nc"
 
! create the grid
 fv3_grid = ESMF_GridCreateMosaic(filename=trim(fname), &
                                  regDecompPTile=decomptile, &
                                  staggerLocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER, &
                                                   ESMF_STAGGERLOC_EDGE1, ESMF_STAGGERLOC_EDGE2/), &
                                  indexflag=ESMF_INDEX_GLOBAL, &
                                  tileFilePath=trim(dir_fix_res), &
                                  rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreateMosaic", ierr)

 end subroutine create_grid_fv3

 subroutine create_grid_gauss(grid_setup, npets, localpet, gauss_grid)

 ! INTENT IN
 type(grid_setup_type), intent(in) :: grid_setup
 integer, intent(in)   :: npets, localpet
 
 ! INTENT OUT 
 type(esmf_grid)                   :: gauss_grid

 type(esmf_polekind_flag)          :: polekindflag(2)

 type(esmf_field)                  :: lon_fields(1)
 type(esmf_field)                  :: lat_fields(1)
 real(esmf_kind_r8), pointer       :: lon_ptr_field(:,:), lon_ptr_coord(:,:)
 real(esmf_kind_r8), pointer       :: lat_ptr_field(:,:), lat_ptr_coord(:,:)
 character(len=15), dimension(1)   :: lon_variables, lat_variables

 integer :: ierr, i_dim, j_dim

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 ! Gaussian increment files lat/lon name
 lon_variables(1) = 'grid_center_lon'
 lat_variables(1) = 'grid_center_lat'
 !lon_variables(1) = 'longitude'
 !lat_variables(1) = 'latitude'

 i_dim = grid_setup%ires
 j_dim = grid_setup%jres

 if (localpet == 0) print*," creating gauss grid for ", i_dim, j_dim
 gauss_grid = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), &
                                    maxIndex=(/i_dim,j_dim/), &
                                    polekindflag=polekindflag, &
                                    periodicDim=1, &
                                    poleDim=2,  &
                                    coordSys=ESMF_COORDSYS_SPH_DEG, &
                                    regDecomp=(/1,npets/),  &
                                    indexflag=ESMF_INDEX_GLOBAL, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridCreate1PeriDim", ierr)

 ! read lat lon coordinates into fields
 call lonlat_read_into_fields(localpet, i_dim, j_dim,  &
                         trim(grid_setup%fname_coord), trim(grid_setup%dir_coord), &
                         grid_setup, gauss_grid, 1, lon_variables, lat_variables, &
                         lon_fields, lat_fields)

 ! add coordinates to the grid
 call ESMF_GridAddCoord(gauss_grid, &
                        staggerloc=ESMF_STAGGERLOC_CENTER, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridAddCoord", ierr)

 nullify(lon_ptr_coord)
 nullify(lon_ptr_field)

 ! get pointer to lon coord
 !if (localpet == 0) print*," getting GridCoord for long"
 print*," getting GridCoord for long", localpet 
 call ESMF_GridGetCoord(gauss_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=1, &
                        farrayPtr=lon_ptr_coord, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord longitude ", ierr)

 ! fetch lon from field into pointer
 call ESMF_FieldGet(lon_fields(1), &
                    farrayPtr=lon_ptr_field, &
                    rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGet longitude", ierr)
 
 ! set coord pointe to field pointer 
 lon_ptr_coord = lon_ptr_field

 nullify(lat_ptr_coord)
 nullify(lat_ptr_field)
 
 ! get pointer to lat coord. Need bounds?
 call ESMF_GridGetCoord(gauss_grid, &
                        staggerLoc=ESMF_STAGGERLOC_CENTER, &
                        coordDim=2, &
                        farrayPtr=lat_ptr_coord, rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
    call error_handler("IN GridGetCoord latitude", ierr)

 ! fetch lat from field into pointer 
 call ESMF_FieldGet(lat_fields(1), &
                    farrayPtr=lat_ptr_field, &
                    rc=ierr)
 if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldGet latitude", ierr)

 ! set lat coord to fiel values 
 lat_ptr_coord =  lat_ptr_field

 ! assign the corners
 ! call ESMF_GridAddCoord(input_grid, &
 !                      staggerloc=ESMF_STAGGERLOC_CORNER, rc=rc)
 !if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
 !   call error_handler("IN GridAddCoord", rc)

 !print*,"- CALL GridGetCoord FOR INPUT GRID X-COORD."
 !nullify(lon_corner_src_ptr)
 !call ESMF_GridGetCoord(input_grid, &
 !                       staggerLoc=ESMF_STAGGERLOC_CORNER, &
 !                       coordDim=1, &
 !                       farrayPtr=lon_corner_src_ptr, rc=rc)
 !if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
 !   call error_handler("IN GridGetCoord", rc)

 !print*,"- CALL GridGetCoord FOR INPUT GRID Y-COORD."
 !nullify(lat_corner_src_ptr)
 !call ESMF_GridGetCoord(input_grid, &
 !                       staggerLoc=ESMF_STAGGERLOC_CORNER, &
 !                       coordDim=2, &
 !                       computationalLBound=clb, &
 !                       computationalUBound=cub, &
 !                       farrayPtr=lat_corner_src_ptr, rc=rc)
 !if(ESMF_logFoundError(rcToCheck=rc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
 !   call error_handler("IN GridGetCoord", rc)

 !do j = clb(2), cub(2)
 !  do i = clb(1), cub(1)
 !    lon_corner_src_ptr(i,j) = longitude(i,1) - (0.5_esmf_kind_r8*deltalon)
 !    if (lon_corner_src_ptr(i,j) > 360.0_esmf_kind_r8) lon_corner_src_ptr(i,j) = lon_corner_src_ptr(i,j) - 360.0_esmf_kind_r8
 !    if (j == 1) then 
 !      lat_corner_src_ptr(i,j) = 90.0_esmf_kind_r8
 !      cycle
 !    endif
 !    if (j == jp1_input) then
 !      lat_corner_src_ptr(i,j) = -90.0_esmf_kind_r8
 !      cycle
 !    endif
 !    lat_corner_src_ptr(i,j) = 0.5_esmf_kind_r8 * (latitude(i,j-1)+ latitude(i,j))
 !  enddo
 !enddo

 ! TO DO destrot the fields 

 end subroutine create_grid_gauss

end  module grids_IO

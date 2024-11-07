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

 public :: setup_grid, &
           read_into_fields, &
           write_from_fields

 contains

!-----------------------------------
! Create ESMF grid objects, with mask if requested
 subroutine setup_grid(localpet, npets, file_type,  & 
                       dir_mask, fname_mask, mask_type, mod_grid, &
                       ! optional arguments for file_type='fv3*'
                       res_atm, dir_fix, &
                       ! optional arguments for file_type='gau*'
                       i_dim, j_dim )

 implicit none

 ! INTENT IN
 character(7), intent(in)       :: file_type
 character(*), intent(in)       :: fname_mask, dir_mask
 character(*), intent(in)       :: mask_type
 integer, intent(in)            :: localpet, npets
 ! NEEDED FOR FV3 GRID
 integer, intent(in), optional            :: res_atm
 character(*), intent(in), optional       :: dir_fix
 ! NEEDED FOR GAUSS GRID
 integer, intent(in), optional            :: i_dim, j_dim

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

 select case (file_type)
 case ('fv3_rst')
     call create_grid_fv3(res_atm, trim(dir_fix), npets, localpet ,mod_grid)
     mask_variable(1) = 'vtype          '
 case ('gau_inc')
     call create_grid_gauss(i_dim, j_dim, npets, localpet, mod_grid)
     mask_variable(1) = 'soilsnow_mask  '
 case default 
     call error_handler("unknown file_type in setup_grid", 1)
 end select

!--------------------------
! Calcalate and add the mask

 if (mask_type=="soil") then 
 ! mask out ocean and glaciers, using vegetation class

     mask_field(1) = ESMF_FieldCreate(mod_grid, &
                                       typekind=ESMF_TYPEKIND_R8, &
                                       staggerloc=ESMF_STAGGERLOC_CENTER, &
                                       name="input variable for mask", &
                                       rc=ierr)
     if(ESMF_logFoundError(rcToCheck=ierr,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,file=__FILE__)) &
        call error_handler("IN FieldCreate, mask_variable", ierr)

     call read_into_fields(localpet, res_atm, res_atm, trim(fname_mask), 1, mask_variable(1), &
                           dir_mask, mask_field(1))

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
 subroutine read_into_fields(localpet, i_dim, j_dim , fname_read, n_vars, variable_list, & 
                             dir_read, fields) 
 
 implicit none 

 ! INTENT IN
 integer, intent(in)             :: localpet, i_dim, j_dim, n_vars
 character(*), intent(in)        :: fname_read
 character(*), intent(in)        :: dir_read

 character(len=15), dimension(n_vars), intent(in)   :: variable_list
 
 ! INTENT OUT 
 type(esmf_field), dimension(n_vars), intent(inout) :: fields

 ! LOCAL
 integer                         :: tt, id_var, ncid, ierr, v
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

 do tt = 1, n_tiles

      ! read from restart
      if (localpet == 0) then

         write(tchar,'(i1)') tt
         fname = dir_read//"/"//fname_read//".tile"//tchar//".nc"

         ierr=nf90_open(trim(fname),NF90_NOWRITE,ncid)
         call netcdf_err(ierr, 'opening: '//trim(fname) )

         do v =1, n_vars
             print *, 'Reading ', trim(variable_list(v)), ' into field, tile', tchar
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

! write variables from ESMF Fields into netcdf restart-like file 

 subroutine write_from_fields(localpet, i_dim, j_dim , fname_out, n_vars, variable_list, & 
                          dir_out, fields) 

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
 character(len=500)     :: fname, dir_fix_res

 if (localpet == 0) print*," creating fv3 grid for ", res_atm 

! pet distribution
 extra = npets / n_tiles
 do tile = 1, n_tiles
   decomptile(:,tile)=(/1,extra/)
 enddo

 ! mosaic file 
 write(rchar,'(i3.3)') res_atm

 dir_fix_res = dir_fix//"/C"//trim(rchar)//"/"

 fname = trim(dir_fix_res)//"/C"//trim(rchar)// "_mosaic.nc"
 
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

 subroutine create_grid_gauss(i_dim, j_dim, npets, localpet, gauss_grid)

 ! INTENT IN
 integer, intent(in)   :: i_dim, j_dim
 integer, intent(in)   :: npets, localpet
 
 ! INTENT OUT 
 type(esmf_grid)                :: gauss_grid

 type(esmf_polekind_flag)         :: polekindflag(2)
 integer :: ierr

 polekindflag(1:2) = ESMF_POLEKIND_MONOPOLE

 if (localpet == 0) print*," creating fv3 gauss for ", i_dim, j_dim
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

 end subroutine create_grid_gauss
end  module grids_IO

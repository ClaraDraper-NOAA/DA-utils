find_package(ESMF 8.0.0 REQUIRED) 

ecbuild_add_executable( TARGET regridStates.x
                        SOURCES regridStates.F90
                                grids_IO.F90
                                utils.F90)

target_link_libraries(
  regridStates.x
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)


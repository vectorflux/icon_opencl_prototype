#----------------------------------------------------------
# Makefile automatically generated by createMakefiles.pl
#----------------------------------------------------------

VPATH = ../../../src/ocean_math_operators:../../../src/atmo_math_operators:../../../src/shared_math_operators:../../../src/shared:../../src:../../../src/atmo_nh_solver:../../../src/drivers:../../../src/OpenCL:../../../src/setup:

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

%.o: %.f
	$(F77) $(F77FLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

%.obj: %.f90
	$(FC) $(FFLAGS) -c $<

SRCS =	mo_solve_nonhydro.f90 mo_nh_global_variables.f90 mo_vertical_grid.f90 \
	mo_nh_driver.f90 mo_nonhydro_state.f90 mo_intp_rbf_coeffs.f90 \
	mo_reshape_arrays.f90 mo_ext_data.f90 mo_intp_coeffs_lsq_bln.f90 \
	mo_interpolation_init.f90 mo_model_domain_import.f90 \
	mo_setup_parameters.f90 mo_domain_geometry_methods.f90 \
	mo_model_domain.f90 mo_global_variables.f90 mo_io_units.f90 \
	mo_impl_constants_grf.f90 mo_datetime.f90 mo_get_utc_date_tr.f90 \
	mo_impl_constants.f90 mo_parallel_ctl.f90 mo_exception.f90 \
	mo_profile.f90 mo_util_string.f90 mo_communication.f90 \
	mo_physical_constants.f90 mo_loopindices.f90 \
	mo_read_netcdf_parallel.f90 mo_namelist.f90 mo_sync.f90 mo_timer.f90 \
	mo_real_timer.f90 mo_math_utilities.f90 mo_mpi.f90 \
	mo_math_constants.f90 mo_kind.f90 mo_hyb_params.f90 \
	mo_atmo_math_operators.f90 mo_atmo_math_gradients.f90 \
	mo_atmo_math_laplace.f90 mo_atmo_math_divrot.f90 \
	mo_oce_math_divrot.f90 mo_interpolation.f90 mo_intp_methods.f90 \
	mo_intp_data_strc.f90 mo_intp_rbf.f90

OBJS =	mo_solve_nonhydro.o mo_nh_global_variables.o mo_vertical_grid.o \
	mo_nh_driver.o mo_nonhydro_state.o mo_intp_rbf_coeffs.o \
	mo_reshape_arrays.o mo_ext_data.o mo_intp_coeffs_lsq_bln.o \
	mo_interpolation_init.o mo_model_domain_import.o \
	mo_setup_parameters.o mo_domain_geometry_methods.o mo_model_domain.o \
	mo_global_variables.o mo_io_units.o mo_impl_constants_grf.o \
	mo_datetime.o mo_get_utc_date_tr.o mo_impl_constants.o \
	mo_parallel_ctl.o mo_exception.o mo_profile.o mo_util_string.o \
	mo_communication.o mo_physical_constants.o mo_loopindices.o \
	mo_read_netcdf_parallel.o mo_namelist.o mo_sync.o mo_timer.o \
	mo_real_timer.o mo_math_utilities.o mo_mpi.o mo_math_constants.o \
	mo_kind.o mo_hyb_params.o mo_atmo_math_operators.o \
	mo_atmo_math_gradients.o mo_atmo_math_laplace.o mo_atmo_math_divrot.o \
	mo_oce_math_divrot.o mo_interpolation.o mo_intp_methods.o \
	mo_intp_data_strc.o mo_intp_rbf.o wrapper.o

all: $(OBJS) icon_testbed

icon_testbed.o: icon_testbed.f90
icon_testbed: icon_testbed.o $(OBJS)
	$(FC) $(LDFLAGS) -o ../bin/$@ $< $(OBJS) $(LIBS)
icon_testbed.exe: 
	$(FC) $(LDFLAGS) -o ../bin/icon_testbed icon_testbed.o $< $(OBJS) $(LIBS)

clean:
	rm -f *.o ../module/*.mod

mo_math_constants.o: mo_math_constants.f90 mo_kind.o

mo_model_domain.o: mo_model_domain.f90 mo_kind.o mo_math_utilities.o \
	mo_communication.o mo_intp_data_strc.o

mo_setup_parameters.o: mo_setup_parameters.f90 mo_global_variables.o \
	mo_namelist.o

mo_exception.o: mo_exception.f90 mo_io_units.o mo_mpi.o mo_kind.o 

mo_intp_coeffs_lsq_bln.o: mo_intp_coeffs_lsq_bln.f90 mo_kind.o \
	mo_math_constants.o mo_physical_constants.o mo_exception.o \
	mo_impl_constants.o mo_model_domain.o mo_math_utilities.o \
	mo_global_variables.o mo_loopindices.o mo_intp_data_strc.o

mo_read_netcdf_parallel.o: mo_read_netcdf_parallel.f90 mo_kind.o mo_mpi.o \
	mo_parallel_ctl.o

mo_intp_rbf.o: mo_intp_rbf.f90 mo_kind.o mo_impl_constants.o \
	mo_model_domain.o mo_global_variables.o mo_loopindices.o \
	mo_intp_data_strc.o

mo_timer.o: mo_timer.f90 mo_real_timer.o mo_global_variables.o

mo_solve_nonhydro.o: mo_solve_nonhydro.f90 mo_kind.o mo_nh_global_variables.o \
	mo_model_domain.o mo_interpolation.o mo_nonhydro_state.o \
	mo_physical_constants.o mo_atmo_math_operators.o mo_vertical_grid.o \
	mo_loopindices.o mo_impl_constants.o mo_impl_constants_grf.o mo_mpi.o \
	mo_intp_data_strc.o mo_global_variables.o mo_impl_constants.o

mo_get_utc_date_tr.o: mo_get_utc_date_tr.f90 mo_kind.o mo_global_variables.o

mo_hyb_params.o: mo_hyb_params.f90 mo_kind.o mo_io_units.o mo_exception.o \
	mo_impl_constants.o mo_global_variables.o mo_physical_constants.o

mo_intp_data_strc.o: mo_intp_data_strc.f90 mo_kind.o mo_impl_constants.o \
	mo_math_utilities.o

mo_atmo_math_laplace.o: mo_atmo_math_laplace.f90 mo_kind.o \
	mo_impl_constants.o mo_intp_methods.o mo_model_domain.o \
	mo_global_variables.o mo_exception.o mo_loopindices.o \
	mo_atmo_math_gradients.o mo_atmo_math_divrot.o

mo_global_variables.o: mo_global_variables.f90 mo_kind.o mo_exception.o \
	mo_impl_constants.o mo_physical_constants.o mo_io_units.o mo_mpi.o \
	mo_datetime.o

mo_model_domain_import.o: mo_model_domain_import.f90 mo_kind.o mo_io_units.o \
	mo_impl_constants.o mo_exception.o mo_model_domain.o \
	mo_physical_constants.o mo_global_variables.o mo_math_utilities.o \
	mo_math_constants.o mo_namelist.o mo_loopindices.o mo_mpi.o \
	mo_io_units.o mo_reshape_arrays.o mo_domain_geometry_methods.o \
	mo_hyb_params.o mo_interpolation_init.o mo_ext_data.o

mo_oce_math_divrot.o: mo_oce_math_divrot.f90 mo_kind.o mo_impl_constants.o \
	mo_interpolation.o mo_model_domain.o mo_global_variables.o \
	mo_exception.o mo_loopindices.o

mo_nh_driver.o: mo_nh_driver.f90 mo_kind.o mo_namelist.o mo_io_units.o \
	mo_nonhydro_state.o mo_nh_global_variables.o mo_ext_data.o \
	mo_model_domain.o mo_interpolation.o mo_interpolation.o mo_datetime.o \
	mo_timer.o mo_exception.o mo_impl_constants.o mo_physical_constants.o \
	mo_math_constants.o mo_math_utilities.o mo_atmo_math_operators.o \
	mo_solve_nonhydro.o mo_mpi.o mo_mpi.o mo_vertical_grid.o \
	mo_loopindices.o

mo_communication.o: mo_communication.f90 mo_kind.o mo_exception.o \
	mo_global_variables.o mo_mpi.o mo_parallel_ctl.o  

mo_loopindices.o: mo_loopindices.f90 mo_kind.o mo_model_domain.o \
	mo_impl_constants.o mo_global_variables.o

mo_reshape_arrays.o: mo_reshape_arrays.f90 mo_kind.o mo_impl_constants.o \
	mo_exception.o mo_model_domain.o mo_physical_constants.o \
	mo_global_variables.o mo_loopindices.o mo_mpi.o

mo_namelist.o: mo_namelist.f90 mo_util_string.o mo_exception.o mo_io_units.o

mo_atmo_math_divrot.o: mo_atmo_math_divrot.f90 mo_kind.o mo_impl_constants.o \
	mo_interpolation.o mo_model_domain.o mo_global_variables.o \
	mo_exception.o mo_loopindices.o

mo_mpi.o: mo_mpi.f90   mo_kind.o mo_io_units.o mo_util_string.o 

mo_physical_constants.o: mo_physical_constants.f90 mo_kind.o

mo_atmo_math_operators.o: mo_atmo_math_operators.f90 mo_atmo_math_gradients.o \
	mo_atmo_math_laplace.o mo_atmo_math_divrot.o

mo_nonhydro_state.o: mo_nonhydro_state.f90 mo_kind.o mo_impl_constants.o \
	mo_exception.o mo_model_domain.o mo_global_variables.o

mo_domain_geometry_methods.o: mo_domain_geometry_methods.f90 mo_kind.o \
	mo_impl_constants.o mo_exception.o mo_model_domain.o \
	mo_physical_constants.o mo_global_variables.o mo_math_utilities.o \
	mo_math_constants.o mo_loopindices.o mo_mpi.o mo_reshape_arrays.o

mo_parallel_ctl.o: mo_parallel_ctl.f90 mo_exception.o mo_io_units.o \
	mo_namelist.o mo_mpi.o mo_mpi.o mo_exception.o

mo_sync.o: mo_sync.f90 mo_kind.o mo_exception.o mo_model_domain.o \
	mo_global_variables.o mo_io_units.o mo_mpi.o mo_parallel_ctl.o \
	mo_communication.o

mo_real_timer.o: mo_real_timer.f90  mo_kind.o mo_exception.o   mo_mpi.o \
	mo_mpi.o mo_parallel_ctl.o mo_kind.o mo_kind.o mo_kind.o

mo_nh_global_variables.o: mo_nh_global_variables.f90 mo_global_variables.o

mo_datetime.o: mo_datetime.f90 mo_kind.o mo_exception.o

mo_vertical_grid.o: mo_vertical_grid.f90 mo_kind.o mo_exception.o \
	mo_model_domain.o mo_ext_data.o mo_global_variables.o \
	mo_setup_parameters.o mo_hyb_params.o mo_impl_constants.o \
	mo_impl_constants_grf.o mo_physical_constants.o \
	mo_atmo_math_operators.o mo_timer.o mo_interpolation.o \
	mo_math_constants.o mo_loopindices.o mo_nonhydro_state.o mo_sync.o \
	mo_parallel_ctl.o

mo_intp_rbf_coeffs.o: mo_intp_rbf_coeffs.f90 mo_kind.o mo_exception.o \
	mo_impl_constants.o mo_model_domain.o mo_math_utilities.o \
	mo_global_variables.o mo_loopindices.o mo_intp_data_strc.o

mo_math_utilities.o: mo_math_utilities.f90 mo_kind.o mo_math_constants.o \
	mo_physical_constants.o mo_exception.o mo_global_variables.o \
	mo_kind.o mo_exception.o mo_kind.o mo_kind.o mo_kind.o

mo_atmo_math_gradients.o: mo_atmo_math_gradients.f90 mo_kind.o \
	mo_impl_constants.o mo_interpolation.o mo_model_domain.o \
	mo_global_variables.o mo_exception.o mo_timer.o mo_loopindices.o

mo_intp_methods.o: mo_intp_methods.f90 mo_kind.o mo_exception.o \
	mo_impl_constants.o mo_model_domain.o mo_global_variables.o \
	mo_loopindices.o mo_intp_data_strc.o mo_timer.o

mo_interpolation.o: mo_interpolation.f90 mo_intp_data_strc.o \
	mo_intp_methods.o mo_intp_rbf.o

mo_interpolation_init.o: mo_interpolation_init.f90 mo_kind.o \
	mo_math_constants.o mo_physical_constants.o mo_impl_constants.o \
	mo_model_domain.o mo_math_utilities.o mo_global_variables.o \
	mo_impl_constants.o mo_intp_data_strc.o mo_loopindices.o \
	mo_exception.o mo_intp_rbf_coeffs.o

icon_testbed.o: icon_testbed.f90 mo_kind.o mo_exception.o mo_namelist.o \
	mo_io_units.o mo_mpi.o mo_timer.o mo_io_units.o mo_datetime.o \
	mo_setup_parameters.o mo_parallel_ctl.o mo_global_variables.o \
	mo_model_domain_import.o mo_profile.o mo_nh_driver.o 

mo_ext_data.o: mo_ext_data.f90 mo_kind.o mo_io_units.o mo_global_variables.o \
	mo_model_domain.o mo_impl_constants.o mo_exception.o \
	mo_reshape_arrays.o mo_intp_methods.o mo_loopindices.o mo_mpi.o \
	mo_parallel_ctl.o mo_communication.o


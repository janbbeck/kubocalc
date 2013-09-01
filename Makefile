# 
# Makefile for kubocalc.x 
# By Lazaro Calderin
# Modded for espresso 4.2
# Jan Beck

PWSCF=../
include $(PWSCF)/make.sys
MODFLAGS       = -I./  -I../Modules  -I../iotk/src \
                 -I../PW/src  -I../PH  -I../EE -I../GIPAW -I../GIPAW/src -I../CPV/src -I../CPV \
                 -I../PW -I../Modules/src -I../PH/src  -I../EE/src
.SUFFIXES :
.SUFFIXES : .o .c .f .f90


MODULES = \
$(PWSCF)/Modules/atom.o \
$(PWSCF)/Modules/autopilot.o \
$(PWSCF)/Modules/basic_algebra_routines.o \
$(PWSCF)/Modules/berry_phase.o \
$(PWSCF)/Modules/bfgs_module.o \
$(PWSCF)/Modules/cell_base.o \
$(PWSCF)/Modules/check_stop.o \
$(PWSCF)/Modules/clocks.o \
$(PWSCF)/Modules/constants.o \
$(PWSCF)/Modules/constraints_module.o \
$(PWSCF)/Modules/control_flags.o \
$(PWSCF)/Modules/descriptors.o \
$(PWSCF)/Modules/dspev_drv.o \
$(PWSCF)/Modules/electrons_base.o \
$(PWSCF)/Modules/error_handler.o \
$(PWSCF)/Modules/environment.o \
$(PWSCF)/Modules/fft_base.o \
$(PWSCF)/Modules/fft_parallel.o \
$(PWSCF)/Modules/fft_scalar.o \
$(PWSCF)/Modules/fft_types.o \
$(PWSCF)/Modules/funct.o \
$(PWSCF)/Modules/input_parameters.o \
$(PWSCF)/Modules/io_files.o \
$(PWSCF)/Modules/io_global.o \
$(PWSCF)/Modules/ions_base.o \
$(PWSCF)/Modules/ions_nose.o \
$(PWSCF)/Modules/kind.o \
$(PWSCF)/Modules/mp_global.o \
$(PWSCF)/Modules/mp_wave.o \
$(PWSCF)/Modules/mp.o \
$(PWSCF)/Modules/mp_base.o \
$(PWSCF)/Modules/metadyn_base.o \
$(PWSCF)/Modules/metadyn_io.o \
$(PWSCF)/Modules/metadyn_vars.o \
$(PWSCF)/Modules/mm_dispersion.o \
$(PWSCF)/Modules/path_formats.o \
$(PWSCF)/Modules/path_variables.o \
$(PWSCF)/Modules/path_opt_routines.o \
$(PWSCF)/Modules/path_io_routines.o \
$(PWSCF)/Modules/path_reparametrisation.o \
$(PWSCF)/Modules/parallel_include.o \
$(PWSCF)/Modules/parameters.o \
$(PWSCF)/Modules/parser.o \
$(PWSCF)/Modules/paw_variables.o \
$(PWSCF)/Modules/printout_base.o \
$(PWSCF)/Modules/pseudo_types.o \
$(PWSCF)/Modules/ptoolkit.o \
$(PWSCF)/Modules/radial_grids.o \
$(PWSCF)/Modules/random_numbers.o \
$(PWSCF)/Modules/read_cards.o \
$(PWSCF)/Modules/read_namelists.o \
$(PWSCF)/Modules/read_ncpp.o \
$(PWSCF)/Modules/read_upf_v1.o \
$(PWSCF)/Modules/read_upf_v2.o \
$(PWSCF)/Modules/read_uspp.o \
$(PWSCF)/Modules/recvec.o \
$(PWSCF)/Modules/splinelib.o \
$(PWSCF)/Modules/stick_base.o \
$(PWSCF)/Modules/task_groups.o \
$(PWSCF)/Modules/timestep.o \
$(PWSCF)/Modules/upf_to_internal.o \
$(PWSCF)/Modules/uspp.o \
$(PWSCF)/Modules/upf.o \
$(PWSCF)/Modules/version.o \
$(PWSCF)/Modules/wannier_new.o \
$(PWSCF)/Modules/wavefunctions.o \
$(PWSCF)/Modules/wave_base.o \
$(PWSCF)/Modules/ws_base.o \
$(PWSCF)/Modules/wrappers.o \
$(PWSCF)/Modules/write_upf_v2.o \
$(PWSCF)/Modules/xml_io_base.o \
$(PWSCF)/Modules/zhpev_drv.o 

EEOBJS = \
$(PWSCF)/EE/ee_mod.o \
$(PWSCF)/EE/ggen_coarse.o \
$(PWSCF)/EE/data_structure_coarse.o \
$(PWSCF)/EE/init_ee.o \
$(PWSCF)/EE/set_fft_dim_coarse.o \
$(PWSCF)/EE/set_mltgrid_dim.o \
$(PWSCF)/EE/write_ee_summary.o \
$(PWSCF)/EE/gcoarse_mod.o

PWOBJS = \
$(PWSCF)/PW/libpw.a 


PPOBJS = \
$(PWSCF)/PP/openfil_pp.o 



all : kubocalc.x

kubocalc.o :
	$(MPIF90) -openmp -nomodule -fpp $(FDFLAGS) $(IFLAGS) $(MODFLAGS) -c $*.f90 -o $*.o

kubocalc.x : kubocalc.o $(PPOBJS) $(EEOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ kubocalc.o $(PPOBJS) $(EEOBJS) $(MODULES) $(PWOBJS) \
	$(LIBOBJS) $(LIBS)
 
clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc

#include make.depend

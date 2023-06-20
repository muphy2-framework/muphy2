
# Possible settings:
# -DSINGLE_PRECISION (single precision if set, double precision if not)
# -DACC_2D_OPT (optimize openacc pragmas for 2D instead of 3D)
SETTINGS =

SRCDIRECTORY = ..

# GCC with OpenMPI or ParaStationMPI
# tested on davinci and JUWELS-booster
# optionally before compiling on davinci use
# module purge; module load gcc mpi/openmpi3-x86_64
# on JUWELS-booster use
# module load GCC ParaStationMPI
    FLAGS = -Ofast -Wall $(SETTINGS) -I$(SRCDIRECTORY)
    CXXFLAGS = $(FLAGS) -std=c++11 -lgfortran -MMD -Wno-unknown-pragmas
    FORTRANFLAGS = $(FLAGS) -ffree-line-length-none
    MPICXX = mpicxx
    MPIF90 = mpifort

# Nvidia/PGI compiler with OpenMPI or ParaStationMPI
# tested on davinci and JUWELS-booster
# before compiling on davinci use
# module purge; module load nvhpc
# on JUWELS-booster use
# module load NVHPC ParaStationMPI
# USE PG-PROFILER: flags:  -Minfo=ccff -Mcuda
# call (openmpi): mpirun <mpi-options> pgprof --cpu-profilng off --openacc-profiling on -o <filename>.%h.%p.%q{OMPI_COMM_WORLD_RANK} ./<executable>
# path to executable in mandatory, './' is sufficient.
    #FLAGS = -fast $(SETTINGS) -I$(SRCDIRECTORY)
    #ACCFLAGS = -acc -ta=tesla,fastmath,keepgpu -Minfo=acc,ccff -Mcuda
    #CXXFLAGS = $(FLAGS) -std=c++11 -pgf90libs -MMD $(ACCFLAGS)
    #FORTRANFLAGS = $(FLAGS) $(ACCFLAGS)
    #MPICXX = mpicxx
    #MPIF90 = mpifort

# Intel compiler with ParaStationMPI
# tested on JUWELS-booster
# before compiling on JUWELS-booster use
# module load Intel ParaStationMPI
    #FLAGS = -Ofast $(SETTINGS) -I$(SRCDIRECTORY)
    #CXXFLAGS = $(FLAGS) -std=c++11 -lifcore -lifport -MMD -Wno-unknown-pragmas
    #FORTRANFLAGS = $(FLAGS)
    #MPICXX = mpicxx
    #MPIF90 = mpifort

# Intel compiler with Intel MPI
# tested on JUWELS-gpu
# before compiling on JUWELS-gpu use
# module load Intel/2019.3.199-GCC-8.3.0 IntelMPI/2019.3.199
    #FLAGS = -Ofast $(SETTINGS) -I$(SRCDIRECTORY)
    #CXXFLAGS = $(FLAGS) -std=c++11 -lifcore -lifport -MMD -Wno-unknown-pragmas
    #FORTRANFLAGS = $(FLAGS)
    #MPICXX = mpiicpc
    #MPIF90 = mpiifort


all: muphy2

XXD:
		@xxd -i ../framework/initial_conditions.cc initial_conditions_files_content.h
		@xxd -i ../framework/parameter.h >> initial_conditions_files_content.h
		@xxd -i ../physics/plasma/setup/setup.F90 >> initial_conditions_files_content.h

OBJS = definitions_f.o utilities.o utilities_f.o initial_conditions.o parameter.o\
		mpi_boundary_f.o mpi_boundary.o restart.o interpolate_f.o interpolate.o\
		converters_physical_quantities_f.o converters_setup_f.o converters_output_f.o converters_models_f.o\
		converters.o setup_f.o setup.o criterion_position.o criterion_j.o criterion_j_ue.o\
		select_scheme.o load_balancing.o\
		output_write.o output.o output_f.o vlasov_pfc_f.o vlasov_lagrange5_f.o vlasov.o fluid_5mom_cweno_f.o\
		fluid_mhd_cweno_f.o mhd.o fluid5.o fluid_10mom_cweno_f.o fluid10.o\
		ohm_f.o ohm.o poisson_f.o poisson.o maxwell_fdtd_f.o maxwell.o VeViM.o VeViMdeltaf.o F10eViM.o F10eViMdeltaf.o\
		MHD1TemperatureOhm.o F10eF10iM.o F5eF10iM.o F5eF5iM.o VeViP.o VeViPdeltaf.o block.o main.o

DEP = $(OBJS:.o=.d)

$(OBJS): | XXD

muphy2: $(OBJS)
	$(MPICXX) -o $@ $(OBJS) $(CXXFLAGS)

-include ${DEP}

definitions_f.o: ../framework/definitions.F90
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

definitions.mod: ../framework/definitions.F90 definitions_f.o
	@true

utilities.o: ../framework/utilities.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

utilities_f.o: ../framework/utilities.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

initial_conditions.o: ../framework/initial_conditions.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

parameter.o: ../framework/parameter.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

mpi_boundary_f.o: ../framework/mpi_boundary.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

mpi_boundary.o: ../framework/mpi_boundary.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

interpolate_f.o: ../framework/interpolate.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

interpolate.o: ../framework/interpolate.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

restart.o: ../framework/restart.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

converters_f.o: ../physics/plasma/converters/converters.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

converters_physical_quantities_f.o: ../physics/plasma/converters/converters_physical_quantities.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

converters_setup_f.o: ../physics/plasma/converters/converters_setup.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

converters_output_f.o: ../physics/plasma/converters/converters_output.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

converters_models_f.o: ../physics/plasma/converters/converters_models.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

converters.o: ../physics/plasma/converters/converters.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

criterion_position.o: ../physics/plasma/criteria/criterion_position.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

criterion_j.o: ../physics/plasma/criteria/criterion_j.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

criterion_j_ue.o: ../physics/plasma/criteria/criterion_j_ue.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

select_scheme.o: ../physics/plasma/criteria/select_scheme.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

load_balancing.o: ../physics/plasma/criteria/load_balancing.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

setup_f.o: ../physics/plasma/setup/setup.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

setup.o: ../physics/plasma/setup/setup.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

output_write.o: ../framework/output_write.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

output_f.o: ../physics/plasma/output/output.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

output.o: ../physics/plasma/output/output.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

vlasov_pfc_f.o: ../physics/plasma/models/vlasov_pfc.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

vlasov_lagrange5_f.o: ../physics/plasma/models/vlasov_lagrange5.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

vlasov.o: ../physics/plasma/models/vlasov.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

fluid_mhd_cweno_f.o: ../physics/plasma/models/fluid_mhd_cweno.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

mhd.o: ../physics/plasma/models/mhd.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

fluid_5mom_cweno_f.o: ../physics/plasma/models/fluid_5mom_cweno.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

fluid5.o: ../physics/plasma/models/fluid5.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

fluid_10mom_cweno_f.o: ../physics/plasma/models/fluid_10mom_cweno.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

fluid10.o: ../physics/plasma/models/fluid10.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

ohm_f.o: ../physics/plasma/models/ohm.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

ohm.o: ../physics/plasma/models/ohm.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

poisson_f.o: ../physics/plasma/models/poisson.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

poisson.o: ../physics/plasma/models/poisson.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

maxwell_fdtd_f.o: ../physics/plasma/models/maxwell_fdtd.F90 definitions.mod
	$(MPIF90) -o $@ -c $< $(FORTRANFLAGS)

maxwell.o: ../physics/plasma/models/maxwell.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

VeViM.o: ../physics/plasma/schemes/VeViM.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

VeViMdeltaf.o: ../physics/plasma/schemes/VeViMdeltaf.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

F10eViM.o: ../physics/plasma/schemes/F10eViM.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

F10eViMdeltaf.o: ../physics/plasma/schemes/F10eViMdeltaf.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

F10eF10iM.o: ../physics/plasma/schemes/F10eF10iM.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

F5eF10iM.o: ../physics/plasma/schemes/F5eF10iM.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

F5eF5iM.o: ../physics/plasma/schemes/F5eF5iM.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

MHD1TemperatureOhm.o: ../physics/plasma/schemes/MHD1TemperatureOhm.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

VeViP.o: ../physics/plasma/schemes/VeViP.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

VeViPdeltaf.o: ../physics/plasma/schemes/VeViPdeltaf.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

block.o: ../framework/block.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

main.o: ../framework/main.cc
	$(MPICXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f *.o *.mod *.d *.gpu initial_conditions_files_content.h

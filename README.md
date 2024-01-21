

# The muphyII Plasma Simulation Framework

![muphyII_reconnection_coupled_simulation_snapshot](images/muphyII_reconnection_coupled_simulation_snapshot.png)

The `muphyII` code is an HPC simulation framework for plasma
physics and targeted at collisionless space plasmas in particular. `muphyII`
utilizes a hierarchy of models with different inherent scales to
address the separation of scales problem typically encountered in
these plasmas. It allows stand-alone use of models as well as a
model-based dynamic and adaptive domain decomposition and resorts to
novel techniques for energy conservation, careful treatment of
inner-domain model boundaries for interface coupling, and robust time
stepping algorithms.


### LICENCE and ACKNOWLEDEMENT

The `muphyII` code is free and open source software licensed under the Mozilla Public
License v. 2.0 which is available at https://mozilla.org/MPL/2.0/.

If the software contributes to findings you decide to present or publish, please be so kind
and cite the reference below. Thank you!

> Allmann-Rahn, F., Lautenbach, S., Deisenhofer, M., Grauer, R., 2023.
*The muphyII code: Multiphysics plasma simulation on large HPC systems*. 
Computer Physics Communications. https://doi.org/10.1016/j.cpc.2023.109064


### MANUAL
User information and documentation is available in the `doc/` folder.


### INSTALLATION

**PRE-COMPILATION**

Make sure to load all required libraries. We use the
Nvidia HPC compiler for GPU compilation and GCC or Intel compiler
for CPU compilation.

Examples are provided in the makefile (`bin/Makefile`).
There are also options to switch from GPU to CPU compilation.
Use the SETTINGS variable in the makefile to switch to
single precision and/or 2D loop optimization.

Uncomment whatever you need or add your own library settings.


**COMPILATION**

Change to bin/ directory and type

`make -j 8`

to compile with 8 processes. If successful, an executable named `muphy2` is generated in the `bin/` folder.

Remove all build files with

`make clean`.


### USE

**SETUP PREPARATION**

Adjust framework/initial_conditions.cc, most importantly:

+ output_directory variable
+ nproc variable, so that nproc[0]*nproc[1]*nproc[2]
is identical to the number of processes you start with
mpirun.

Other setups are available in the
`physics/plasma/setup/initial_conditions/` directory --
simply copy the respective file to `framework/initial_conditions.cc`.

Changes to the parameter file require recompilation.

**Employment**

Run the `muphy2` executive in parallel using mpirun, srun, or your sumbission manager. Examples are 

`mpirun -np 16 --oversubscribe ./muphy2`

`srun ./muphy2`


### USER SUPPORT

Contributions, suggestions, and notes on issues are always welcome. Unfortunately, we do
not have the ressources to offer individual user support.


### FILE STRUCTURE OVERVIEW

+ `bin/`: Contains the Makefile and after compilation the executable and build files.
+ `framework/`: Files that are not limited to a special
use case (such as plasma simulations), e.g.: main, parameter, block, MPI/boundary, file writing, interpolation, ...
+ `physics/*/models/`: Physical models, i.e. their data, and functions implementing the respective numerical solvers.
+ `physics/*/schemes/`: Schemes that determine the execution order of the models' solver functions.
+ `physics/*/criteria/`: Criteria that determine which scheme is run in which block and functions for
allocating schemes.
+ `physics/*/converters/`: Functions that convert physical quantities or models into each other, or convert
initial conditions to models and models to output.
+ `physics/*/output/`: Preparation of output and call of file writing functions.
+ `physics/*/setup/`: Initial configurations. In setup.F90, typical setups are defined that
can be further specified in the framework/initial_conditions.cc file.


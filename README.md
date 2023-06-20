# The muphyII Plasma Simulation Framework

### PURPOSE

The `muphyII` code is a simulation framework for collisionless plasma
physics and targeted at space plasmas in particular. The `muphyii`
code utilizes a hierarchy of models with different inherent scales to
address the separation of scales problem typically encountered in
these plasmas. It allows stand-alone use of models as well as a
model-based dynamic and adaptive domain decomposition and resorts to
novel techniques for energy conservation, careful treatment of
inner-domain model boundaries for interface coupling, and robust time
stepping algorithms.


### LICENCE and ACKNOWLEDEMENT

The `muphyII` code is published under MPL-2.0

If this software contributes to findings you decide to present or
publish, we appreciate you citing the reference below. Thank you!

> Allmann-Rahn, F., Lautenbach, S., Deisenhofer, M., Grauer, R., 2023. *The muphyII Code: Multiphysics Plasma Simulation on Large
HPC systems*. arXiv. https://doi.org/10.48550/arXiv.2305.01487


### INSTALLATION

**PRE-COMPILATION**

Make sure to load all required libraries. It is recommended to use the
Nvidia HPC compiler with CUDA and a suitable MPI-3.0 implementation.

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

`mpirun -np 8 --machinefiles/machines-c2-gpu ./muphy2`

`mpirun -np 128 --machinefiles/machines-all ./muphy2`

`srun ./muphy2`


### USER SUPPORT

Contributions, suggestions, and notes on issues are always welcome. Unfortunately, we do
not have the ressources to offer individual user support.


### FRAMEWORK 101

We recommend every prospective user read the user guide located in the `doc/` folder.

**FILE STRUCTURE**

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

Currently, each block corresponds to one MPI process and holds one scheme.
Each scheme holds one or multiple models.


**LANGUAGES**

C++ is used for the program logic, Fortran for calculations.
For C++, we loosely follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
For Fortran, the [Fortran90 Best Practices](https://www.fortran90.org/src/best-practices.html) is
a helpful resource. Don't forget to put the _PRC appendix on Fortran floating point numbers
to get the correct precision.


**TESTING**

So far no explicit test files are available -- one has to keep track by oneself which
code parts could potentially be affected by changes and test them manually.
Typically, a good test setup is the one in physics/plasma/setup/parameters/parameter_gem.cc
which runs several schemes and the coupling between them.



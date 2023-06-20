/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef PARAMETER_H_
#define PARAMETER_H_

#include <string>
#include <unordered_map>
#include "framework/definitions.h"
#include "framework/criterion.h"

struct Species  // species in the plasma (electrons, ions, ...)
{
    real m,                     // mass
        q,                      // charge
        T0,                     // initial temperature
        vb[3],                  // beginning of velocity space
        ve[3],                  // end of velocity space
        dv[3],                  // grid spacing velocity space
        ten_moment_k0;          // relaxation parameter in closure to ten moment equations
    int res_v[3];               // velocity space resolution
    long int ncells_xv_nobd,    // number of cells in phase space without boundary
        ncells_xv;              // number of cells in phase space with boundary
};

struct Parameter {

    Parameter();
    ~Parameter();
    void init();
    void write_initial_conditions_files();

    // program flow
    bool restart = false;
    int nupdates_scheme = 0;
    std::string setup;
    std::string output_directory;
    int noutputs_vtk = 100, noutputs_csv = 100;
    bool restart_output = true;
    std::string phase_space_output = "none"; // may be: "none", "slice", "full", "quarter"
    int phase_space_output_direction = 0; // for phase_space_output == "slice"
    // if for example phase_space_output_direction == 0 (x direction), then
    // phase_space_output_slice_coordinates sets the indices y, z, v_y, v_z
    // of the distribution function, where the x-v_x slice will be taken
    int phase_space_output_slice_index[4] = {0, 0, 0, 0};
    bool heatflux_output = false;
    int nconsecutive_outputs = 1; // if set to e.g. 3, then 3 time steps in a row are output for calculation of d/dt
                                  // better not be combined with electron_subcycling==true
    bool freeze_dt = false;

    // domain and resolution
    int nproc[3],               // number of processes in each direction
        res_x_total[3],         // position space resolution (number of cells) in each direction
        res_x[3],               // number of cells per block, without boundary
        res_x_minus_one[3],     // for fortran code compatibility
        bd[3],                  // number of boundary cells in each direction
        dimensionality_x,       // number of dimensions in position space
        dimensionality_v,       // number of dimensions in velocity space
        ncells_x_nobd,          // number of cells in position space without boundary
        ncells_x,               // number of cells in position space with boundary
        nslots_per_node = -1;   // number of CPU cores or GPUs per node; this needs to be set for load balancing
                                // when there are more MPI processes than CPU cores or GPUs (using e.g. CUDA MPS)
    real xb[3],                 // beginning of position space
        xe[3],                  // end of position space
        dx[3],                  // grid spacing
        xb_loc[3] = {0.,0.,0.}; // beginning of block in position space

    // time
    real t_begin = 0.,          // time begin
        t_end;                  // time end

    // physical constants and models
    int default_scheme = -1;
    enum model_names {kMaxwell, kPoisson, kOhm, kVlasov, kFluid10, kFluid5, kMHD1Temperature};
    enum scheme_names {kVeViM, kF10eViM, kF10eF10iM, kF5eF10iM, kF5eF5iM, kMHD1TemperatureOhm,
                       kVeViP, kVeViPdeltaf, kVeViMdeltaf, kF10eViMdeltaf};
    int nschemes = 5;
    // enum integer values in this array must be ordered from low to high
    int scheme_hierarchy[5] = {kVeViMdeltaf, kF10eViMdeltaf, kF10eF10iM, kF5eF10iM, kF5eF5iM};
    criterion::Criterion* criterion = nullptr;
    enum criterion_names {kCriterionPosition, kCriterionJ, kCriterionJUe};
    std::string vlasov_solver = "pfc";  // may be: "pfc", "lagrange5"
    std::string poisson_solver = "sor"; // may be: "sor", "gauss_seidel"
    // x lower, x upper, y lower, y upper, z lower, z upper
    enum direction_names {kXL, kXU, kYL, kYU, kZL, kZU, kNDirections};
    real c0, mu0, eps0;                             // speed of light, vacuum permeability, vacuum permittivity
    int nspecies;                                   // number of species in the plasma
    enum species_names {kElectron,kIon};            // available species
    std::string species_names_str[2] = {"e","i"};   // available species, string for output
    Species* species;

    // setup specific
    std::unordered_map<std::string,real> setup_var;
    std::unordered_map<std::string,bool> setup_var_bool;
    std::unordered_map<std::string,int> setup_var_int;

    // numerics
    real cfl,                           // CFL condition cfl=v_max*dt/dx
      cfl_maxwell = 0.95,               // of FDTD maxwell solver
      cfl_vlasov_ions = -1.,            // of Vlasov solver ions; not used if < 0
                                        // especially useful for hybrid fluid electrons, Vlasov ions schemes
      cfl_fluid_ions = -1.,             // of fluid solver ions; not used if < 0
                                        // especially useful for hybrid fluid electrons, Vlasov-deltaf ions scheme
      cweno_limiter = 1.e-5,            // minimum allowed density/temperature
      cweno_epsilon = 1.e-5,
      cweno_epsilon_deltaf = -1.,       // separate cweno_epsilon for deltaf schemes (useful for F10eViMdeltaf scheme
                                        // to use smaller cweno_epsilon for the Vlasov ions; cannot be used in coupling);
                                        // if not set (i.e. < 0), will be set to be equal to cweno_epsilon
      mhd_resistivity = 0.,
      poisson_schwarz_iterations = 1000,
      poisson_convergence_threshold = 1.e-7;

    std::string ten_moment_closure;     // may be: "gradient_P", "gradient_T", "isotropization", "gradient_sym", "none"
    int initial_maxwell_steps_per_fluid_step = 21,
        gradient_closure_subcycles = 8;
    bool relax_fdtd_cfl = false,                    // if true, maxwell FDTD cfl condition only has to be
                                                    // fulfilled for maxwell subcycling dt, not for plasma dt
         cweno_calc_j_from_flux = true,             // better fulfillment of Gauss's law, but not applicable in pure Vlasov models
         ten_moment_vlasov_coupling_fit = true,     // improves spatial coupling of ten moment and Vlasov models, but
                                                    // increases computational cost and may reduce stability
         five_ten_moment_coupling_fit = true,       // improves spatial coupling of five and ten moment models, but
                                                    // may reduce stability of e.g. the gradient_P closure
         electron_subcycling = false;               // if true, advance ions only on ion time scales


    // boundary conditions
    bool is_periodic[3] = {true, true, true};
    // order: x low, x high, y low, y high, z low, z high
    // p periodic, q antiperiodic, m mirrored, s symmetric, a antisymmetric, z zero
    // example: 'ppmmpa'
    std::string bd_cond_v[3],           // velocity in each direction
            bd_cond_E[3],
            bd_cond_B[3],
            bd_cond_phi[1],             // electric potential
            bd_cond_n[1],
            bd_cond_nuu[1],             // n u \dot u
            bd_cond_P[6],               // symmetric pressure tensor P_xx, P_xy, P_xz, P_yy, P_yz, P_zz
            bd_cond_Q_h[10],            // symmetric heatflux tensor Q_xxx, Q_xxy, ..., Q_zzz
            bd_cond_five_moments[5],
            bd_cond_ten_moments[10];
};

#endif // PARAMETER_H_

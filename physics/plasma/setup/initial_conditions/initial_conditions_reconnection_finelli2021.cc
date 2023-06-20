/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "parameter.h"
#include <math.h>
#include <fstream>

#include "bin/initial_conditions_files_content.h"


void Parameter::init() {

  // Single current sheet version of the setup in
  // https://doi.org/10.1051/0004-6361/202140279
  // Finelli et al. "Bridging hybrid- and full-kinetic models with Landau-fluid electrons"

  default_scheme = kF10eF10iM;
  //default_scheme = kF10eViMdeltaf;

  output_directory = "./output/";
  noutputs_vtk = 75;
  noutputs_csv = 75;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 1./100.;
  species[kIon].m = 1.;
  species[kElectron].T0 = 1. / 10.;
  species[kIon].T0 = 4. / 10.;
  c0 = 20.;
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "harris_sheet";
  setup_var["n_bg"] = 1.;
  setup_var["T_bg_e"] = species[kElectron].T0;
  setup_var["T_bg_i"] = species[kIon].T0;
  setup_var_bool["drifting_background"] = false;
  setup_var_bool["sine_perturbation"] = false;
  setup_var["lambda"] = 0.85;
  setup_var["psi"] = 0.02;
  setup_var["guide_field"] = .25;
  setup_var["noise_level"] = 0.;

  // domain
  nproc[0] = 4;
  nproc[1] = 2;
  nproc[2] = 1;

  res_x_total[0] = 1024;
  res_x_total[1] = 512;
  res_x_total[2] = 1;

  species[kElectron].res_v[0] = 36;
  species[kElectron].res_v[1] = 36;
  species[kElectron].res_v[2] = 36;
  species[kIon].res_v[0] = 16;
  species[kIon].res_v[1] = 16;
  species[kIon].res_v[2] = 16;

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 2;
  dimensionality_v = 3;

  xb[0] = -12.*M_PI;
  xb[1] = -6.*M_PI;
  xb[2] = -1.;

  xe[0] = 12.*M_PI;
  xe[1] = 6.*M_PI;
  xe[2] = 1.;

  species[kElectron].vb[0] = -20.;
  species[kElectron].vb[1] = -20.;
  species[kElectron].vb[2] = -20.;

  species[kElectron].ve[0] = 20.;
  species[kElectron].ve[1] = 20.;
  species[kElectron].ve[2] = 20.;

  species[kIon].vb[0] = -4.;
  species[kIon].vb[1] = -4.;
  species[kIon].vb[2] = -4.;

  species[kIon].ve[0] = 4.;
  species[kIon].ve[1] = 4.;
  species[kIon].ve[2] = 4.;

  // time
  t_end = 150.;

  // numerics
  cfl = 0.1;
  cfl_fluid_ions = 0.2;
  cfl_vlasov_ions = 0.5;
  ten_moment_closure = "gradient_T";
  species[kElectron].ten_moment_k0 = 1./(2.*sqrt(2./M_PI)) / sqrt(species[kElectron].m);
  species[kIon].ten_moment_k0 = 1./(2.*sqrt(2./M_PI)) / sqrt(species[kIon].m);
  gradient_closure_subcycles = 8; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 7;
  electron_subcycling = true;
  cweno_epsilon = 0.1;

  // boundary conditions
  is_periodic[0] = true;
  is_periodic[1] = false;
  is_periodic[2] = true;

  bd_cond_v[0] = "ppsspp";
  bd_cond_v[1] = "ppzzpp";
  bd_cond_v[2] = "ppsspp";

  bd_cond_E[0] = "ppzzpp";
  bd_cond_E[1] = "ppsspp";
  bd_cond_E[2] = "ppzzpp";

  bd_cond_B[0] = "ppsspp";
  bd_cond_B[1] = "ppzzpp";
  bd_cond_B[2] = "ppsspp";

  bd_cond_n[0] = "ppsspp";
  bd_cond_nuu[0] = "ppsspp";

  bd_cond_P[0] = "ppsspp"; // xx
  bd_cond_P[1] = "ppzzpp"; // xy
  bd_cond_P[2] = "ppsspp"; // xz
  bd_cond_P[3] = "ppsspp"; // yy
  bd_cond_P[4] = "ppzzpp"; // yz
  bd_cond_P[5] = "ppsspp"; // zz

  bd_cond_Q_h[0] = "ppsspp"; // xxx
  bd_cond_Q_h[1] = "ppzzpp"; // xxy
  bd_cond_Q_h[2] = "ppsspp"; // xxz
  bd_cond_Q_h[3] = "ppsspp"; // xyy
  bd_cond_Q_h[4] = "ppzzpp"; // xyz
  bd_cond_Q_h[5] = "ppsspp"; // xzz
  bd_cond_Q_h[6] = "ppzzpp"; // yyy
  bd_cond_Q_h[7] = "ppsspp"; // yyz
  bd_cond_Q_h[8] = "ppzzpp"; // yzz
  bd_cond_Q_h[9] = "ppsspp"; // zzz
}


void Parameter::write_initial_conditions_files() {

    std::ofstream os((output_directory + "/log/initial_conditions_cc.txt").c_str(), std::ios::out);
    os<<___framework_initial_conditions_cc;
    os.close();
    os.clear();
    os.open((output_directory + "/log/parameter_h.txt").c_str(), std::ios::out);
    os<<___framework_parameter_h;
    os.close();
    os.clear();
    os.open((output_directory + "/log/setup_F90.txt").c_str(), std::ios::out);
    os<<___physics_plasma_setup_setup_F90;
    os.close();
}


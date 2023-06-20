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

#include "physics/plasma/criteria/criterion_j.h"


void Parameter::init() {

  // GEM

  default_scheme = kF5eF5iM;
  real thresholds[4] = {1.25, 0.75, 0.35, 0.2};
  criterion = new criterion::Criterion_j(nschemes, scheme_hierarchy, thresholds);
  nupdates_scheme = 40;

  output_directory = "./output/";
  noutputs_vtk = 40;
  noutputs_csv = 40;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 0.04;
  species[kIon].m = 1.;
  species[kElectron].T0 = 1. / 12.;
  species[kIon].T0 = 5. / 12.;
  c0 = 20.;
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "harris_sheet";
  setup_var["n_bg"] = 0.2;
  setup_var["T_bg_e"] = species[kElectron].T0;
  setup_var["T_bg_i"] = species[kIon].T0;
  setup_var_bool["drifting_background"] = true;
  setup_var_bool["sine_perturbation"] = true;
  setup_var["lambda"] = 0.5;
  setup_var["psi"] = 0.1;
  setup_var["guide_field"] = 0.;
  setup_var["noise_level"] = 0.;

  // domain
  nproc[0] = 8;
  nproc[1] = 16;
  nproc[2] = 1;

  res_x_total[0] = 128;
  res_x_total[1] = 64;
  res_x_total[2] = 1;

  for (int i = 0; i < nspecies; ++i) {
    species[i].res_v[0] = 16; // better 32
    species[i].res_v[1] = 16; // better 32
    species[i].res_v[2] = 16; // better 32
  }

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 2;
  dimensionality_v = 3;

  xb[0] = 0.;
  xb[1] = -2. * M_PI;
  xb[2] = -1.;

  xe[0] = 4. * M_PI;
  xe[1] = 0.;
  xe[2] = 1.;

  species[kElectron].vb[0] = -10.;
  species[kElectron].vb[1] = -10.;
  species[kElectron].vb[2] = -10.;

  species[kElectron].ve[0] = 10.;
  species[kElectron].ve[1] = 10.;
  species[kElectron].ve[2] = 10.;

  species[kIon].vb[0] = -5.;
  species[kIon].vb[1] = -5.;
  species[kIon].vb[2] = -5.;

  species[kIon].ve[0] = 5.;
  species[kIon].ve[1] = 5.;
  species[kIon].ve[2] = 5.;

  // time
  t_end = 40.;

  // numerics
  cfl = 0.1;
  ten_moment_closure = "gradient_P";
  //species[kElectron].ten_moment_k0 = 1./3. / sqrt(species[kElectron].m);
  //species[kIon].ten_moment_k0 = 1./3. / sqrt(species[kIon].m);
  // best agreement with Vlasov for small system:
  species[kElectron].ten_moment_k0 = 3.;
  species[kIon].ten_moment_k0 = 3.;
  gradient_closure_subcycles = 4; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 15;

  // boundary conditions
  is_periodic[0] = false;
  is_periodic[1] = false;
  is_periodic[2] = true;

  bd_cond_v[0] = "aasspp";
  bd_cond_v[1] = "sszapp";
  bd_cond_v[2] = "sssspp";

  bd_cond_E[0] = "aazspp";
  bd_cond_E[1] = "sssapp";
  bd_cond_E[2] = "sszspp";
                     
  bd_cond_B[0] = "sssapp";
  bd_cond_B[1] = "aazspp";
  bd_cond_B[2] = "aaaapp";

  bd_cond_n[0] = "sssspp";
  bd_cond_nuu[0] = "sssspp";

  bd_cond_P[0] = "sssspp";
  bd_cond_P[1] = "aazapp";
  bd_cond_P[2] = "aasspp";
  bd_cond_P[3] = "sssspp";
  bd_cond_P[4] = "sszapp";
  bd_cond_P[5] = "sssspp";

  bd_cond_Q_h[0] = "aasspp"; // xxx
  bd_cond_Q_h[1] = "sszapp"; // xxy
  bd_cond_Q_h[2] = "sssspp"; // xxz
  bd_cond_Q_h[3] = "aasspp"; // xyy
  bd_cond_Q_h[4] = "aazapp"; // xyz
  bd_cond_Q_h[5] = "aasspp"; // xzz
  bd_cond_Q_h[6] = "sszapp"; // yyy
  bd_cond_Q_h[7] = "sssspp"; // yyz
  bd_cond_Q_h[8] = "sszapp"; // yzz
  bd_cond_Q_h[9] = "sssspp"; // zzz
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


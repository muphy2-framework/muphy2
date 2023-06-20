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

  // island coalescence

  default_scheme = kF10eF10iM;

  output_directory = "./output/";
  noutputs_vtk = 100;
  noutputs_csv = 100;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 0.04;
  species[kIon].m = 1.;
  species[kElectron].T0 = 0.25;
  species[kIon].T0 = 0.25;
  c0 = 15.;
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "island_coalescence";
  setup_var["lambda"] = 5.;
  setup_var["eta"] = 0.4;
  setup_var["psi"] = 0.1;
  setup_var["n_bg"] = 0.2;
  setup_var["guide_field"] = 0.;
  setup_var["noise_level"] = 1.e-4;

  // domain
  nproc[0] = 2;
  nproc[1] = 4;
  nproc[2] = 1;

  res_x_total[0] = 220;
  res_x_total[1] = 440;
  res_x_total[2] = 44;

  for (int i = 0; i < nspecies; ++i) {
    species[i].res_v[0] = 16;
    species[i].res_v[1] = 16;
    species[i].res_v[2] = 16;
  }

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 3;
  dimensionality_v = 3;

  xb[0] = -setup_var["lambda"]*M_PI;
  xb[1] = -2.*setup_var["lambda"]*M_PI;
  xb[2] = -1.*M_PI;

  xe[0] = setup_var["lambda"]*M_PI;
  xe[1] = 2.*setup_var["lambda"]*M_PI;
  xe[2] = 1.*M_PI;

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
  t_end = xe[1]*5.; // 2.5 t_A = 2.5 * L_y

  // numerics
  cfl = 0.1;
  ten_moment_closure = "gradient_P";
  species[kElectron].ten_moment_k0 = 1./3. / sqrt(species[kElectron].m);
  species[kIon].ten_moment_k0 = 1./3. / sqrt(species[kIon].m);
  gradient_closure_subcycles = 8; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 15;

  // boundary conditions
  is_periodic[0] = false;
  is_periodic[1] = true;
  is_periodic[2] = true;

  bd_cond_v[0] = "zzpppp";
  bd_cond_v[1] = "mmpppp";
  bd_cond_v[2] = "mmpppp";

  bd_cond_E[0] = "zzpppp";
  bd_cond_E[1] = "zzpppp";
  bd_cond_E[2] = "zzpppp";

  bd_cond_B[0] = "zzpppp";
  bd_cond_B[1] = "qqpppp";
  bd_cond_B[2] = "zzpppp"; // adjust if guide field

  bd_cond_n[0]    = "mmpppp";
  bd_cond_nuu[0]  = "mmpppp";

  bd_cond_P[0] = "mmpppp";
  bd_cond_P[1] = "zzpppp";
  bd_cond_P[2] = "zzpppp";
  bd_cond_P[3] = "mmpppp";
  bd_cond_P[4] = "mmpppp";
  bd_cond_P[5] = "mmpppp";

  bd_cond_Q_h[0] = "zzpppp"; // xxx
  bd_cond_Q_h[1] = "mmpppp"; // xxy
  bd_cond_Q_h[2] = "zmpppp"; // xxz
  bd_cond_Q_h[3] = "zzpppp"; // xyy
  bd_cond_Q_h[4] = "mzpppp"; // xyz
  bd_cond_Q_h[5] = "zzpppp"; // xzz
  bd_cond_Q_h[6] = "mmpppp"; // yyy
  bd_cond_Q_h[7] = "zmpppp"; // yyz
  bd_cond_Q_h[8] = "mmpppp"; // yzz
  bd_cond_Q_h[9] = "zmpppp"; // zzz
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


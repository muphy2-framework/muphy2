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

  // reconnection ipic3D setup

  default_scheme = kF10eF10iM;

  output_directory = "./output/";
  noutputs_vtk = 50;
  noutputs_csv = 50;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 1./256.;
  species[kIon].m = 1.;
  species[kElectron].T0 = 1. / 12.;
  species[kIon].T0 = 5. / 12.;
  c0 = 40.;
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "double_harris_sheet_gaussian_hump";
  setup_var["n_bg"] = 0.1;
  setup_var["lambda"] = 0.5;
  setup_var["guide_field"] = 0.1;
  setup_var["noise_level"] = 0.;

  // domain
  nproc[0] = 2;
  nproc[1] = 4;
  nproc[2] = 1;

  res_x_total[0] = 768;
  res_x_total[1] = 1024;
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
  xb[1] = 0.;
  xb[2] = -1.;

  xe[0] = 30.;
  xe[1] = 40.;
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
  t_end = 10.;

  // numerics
  cfl = 0.1;
  cweno_limiter = 0.0001;
  cweno_epsilon = 1.e-5;
  ten_moment_closure = "gradient_T";
  species[kElectron].ten_moment_k0 = 1./(2.*sqrt(2./M_PI)) / sqrt(species[kElectron].m);
  species[kIon].ten_moment_k0 = 1./(2.*sqrt(2./M_PI)) / sqrt(species[kIon].m);
  gradient_closure_subcycles = 8; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 15;

  // boundary conditions
  is_periodic[0] = true;
  is_periodic[1] = true;
  is_periodic[2] = true;

  bd_cond_v[0] = "pppppp";
  bd_cond_v[1] = "pppppp";
  bd_cond_v[2] = "pppppp";

  bd_cond_E[0] = "pppppp";
  bd_cond_E[1] = "pppppp";
  bd_cond_E[2] = "pppppp";

  bd_cond_B[0] = "pppppp";
  bd_cond_B[1] = "pppppp";
  bd_cond_B[2] = "pppppp";

  bd_cond_n[0] = "pppppp";
  bd_cond_nuu[0] = "pppppp";

  bd_cond_P[0] = "pppppp";
  bd_cond_P[1] = "pppppp";
  bd_cond_P[2] = "pppppp";
  bd_cond_P[3] = "pppppp";
  bd_cond_P[4] = "pppppp";
  bd_cond_P[5] = "pppppp";

  bd_cond_Q_h[0] = "pppppp"; // xxx
  bd_cond_Q_h[1] = "pppppp"; // xxy
  bd_cond_Q_h[2] = "pppppp"; // xxz
  bd_cond_Q_h[3] = "pppppp"; // xyy
  bd_cond_Q_h[4] = "pppppp"; // xyz
  bd_cond_Q_h[5] = "pppppp"; // xzz
  bd_cond_Q_h[6] = "pppppp"; // yyy
  bd_cond_Q_h[7] = "pppppp"; // yyz
  bd_cond_Q_h[8] = "pppppp"; // yzz
  bd_cond_Q_h[9] = "pppppp"; // zzz

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


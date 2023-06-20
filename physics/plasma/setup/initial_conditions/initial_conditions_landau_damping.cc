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

  // Landau damping

  default_scheme = kVeViP;

  output_directory = "./output/";
  noutputs_vtk = 100;
  noutputs_csv = 1000;

  phase_space_output = "slice";

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 1.;
  species[kIon].m = 1836.;
  species[kElectron].T0 = 1.;
  species[kIon].T0 = 1.;
  c0 = 1.;
  mu0 = 1./(c0*c0);
  eps0 = 1.;

  // setup specific
  setup = "landau_damping";
  setup_var["alpha"] = 0.01;

  // domain
  nproc[0] = 8;
  nproc[1] = 1;
  nproc[2] = 1;

  res_x_total[0] = 128;
  res_x_total[1] = 1;
  res_x_total[2] = 1;

  for (int i = 0; i < nspecies; ++i) {
    species[i].res_v[0] = 128;
    species[i].res_v[1] = 1;
    species[i].res_v[2] = 1;
  }

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 1;
  dimensionality_v = 1;

  xb[0] = -2.*M_PI;
  xb[1] = -2.*M_PI;
  xb[2] = -2.*M_PI;

  xe[0] = 2.*M_PI;
  xe[1] = 2.*M_PI;
  xe[2] = 2.*M_PI;

  species[kElectron].vb[0] = -6.;
  species[kElectron].vb[1] = -6.;
  species[kElectron].vb[2] = -6.;

  species[kElectron].ve[0] = 6.;
  species[kElectron].ve[1] = 6.;
  species[kElectron].ve[2] = 6.;

  species[kIon].vb[0] = -6.;
  species[kIon].vb[1] = -6.;
  species[kIon].vb[2] = -6.;

  species[kIon].ve[0] = 6.;
  species[kIon].ve[1] = 6.;
  species[kIon].ve[2] = 6.;

  // time
  t_end = 50.;

  // numerics
  cfl = 0.9; // lower for fluid
  vlasov_solver = "pfc";
  ten_moment_closure = "gradient_T";
  species[kElectron].ten_moment_k0 = 1./8.;
  species[kIon].ten_moment_k0 = 1./8.;
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

  bd_cond_phi[0] = "pppppp";

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


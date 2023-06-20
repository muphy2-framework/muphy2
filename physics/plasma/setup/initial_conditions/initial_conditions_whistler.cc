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

#include "physics/plasma/criteria/criterion_position.h"


void Parameter::init() {

  // whistler wave

  output_directory = "./output/";
  noutputs_vtk = 10;
  noutputs_csv = 10;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 1./25.;
  species[kIon].m = 1.;
  species[kElectron].T0 = 5.12e-4;
  species[kIon].T0 = 5.12e-4;
  c0 = 20.;
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "whistler_wave";

  // domain
  nproc[0] = 16;
  nproc[1] = 8;
  nproc[2] = 1;

  res_x_total[0] = 64;
  res_x_total[1] = 32;
  res_x_total[2] = 1;

  for (int i = 0; i < nspecies; ++i) {
    species[i].res_v[0] = 16;
    species[i].res_v[1] = 16;
    species[i].res_v[2] = 16;
  }

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 2;
  dimensionality_v = 3;

  xb[0] = -11.18033989;
  xb[1] = -5.590169945;
  xb[2] = -1.;

  xe[0] = 11.18033989;
  xe[1] = 5.590169945;
  xe[2] = 1.;

  species[kElectron].vb[0] = -0.9;
  species[kElectron].vb[1] = -0.9;
  species[kElectron].vb[2] = -0.9;

  species[kElectron].ve[0] = 0.9;
  species[kElectron].ve[1] = 0.9;
  species[kElectron].ve[2] = 0.9;

  species[kIon].vb[0] = -0.18;
  species[kIon].vb[1] = -0.18;
  species[kIon].vb[2] = -0.18;

  species[kIon].ve[0] = 0.18;
  species[kIon].ve[1] = 0.18;
  species[kIon].ve[2] = 0.18;

  // criterion
  default_scheme = kF5eF5iM;
  int base_point_x[4] = {nproc[0]/2-1, nproc[0]/2-2, nproc[0]/2-3, nproc[0]/2-4},
      base_point_y[4] = {nproc[1]/2-1, nproc[1]/2-2, nproc[1]/2-3, nproc[1]/2-4},
      base_point_z[4] = {0, 0, 0, 0},
      width[4] =        {2, 4, 6, 8},
      height[4] =       {2, 4, 6, 8},
      depth[4] =        {1, 1, 1, 1};
  criterion = new criterion::Criterion_position(nschemes, scheme_hierarchy,
                                                base_point_x, base_point_y, base_point_z,
                                                width, height, depth);

  // time
  t_end = 7.31;

  // numerics
  cfl = 0.005;
  ten_moment_closure = "gradient_P";
  species[kElectron].ten_moment_k0 = 1./3. /sqrt(species[kElectron].m);
  species[kIon].ten_moment_k0 = 1./3. /sqrt(species[kIon].m);
  gradient_closure_subcycles = 8; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 15;
  ten_moment_vlasov_coupling_fit = true;
  five_ten_moment_coupling_fit = true;

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


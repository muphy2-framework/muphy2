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

  // Orszag-Tang vortex

  // ATTENTION: m = 1/25

  default_scheme = kF10eF10iM;

  output_directory = "./output/";
  noutputs_vtk = 50;
  noutputs_csv = 50;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 0.04; // Groselj et al.: 0.01
  species[kIon].m = 1.;
  // A1, A2
  species[kElectron].T0 = .05;
  species[kIon].T0 = .05;
  // B1, B2
  //species[kElectron].T0 = .25;
  //species[kIon].T0 = .25;
  c0 = 18.174; // A1, A2
  //c0 = 38.222; // B1, B2
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "orszag_tang";
  setup_var["n_bg"] = 1.;
  setup_var["L"] = 8. * M_PI; // Groselj et al.: 8. * M_PI
  setup_var["delta_u"] = 0.2; // A1
  //setup_var["delta_u"] = 0.1; // A2
  //setup_var["delta_u"] = 0.3; // B1
  //setup_var["delta_u"] = 0.15; // B2

  // domain
  nproc[0] = 16;
  nproc[1] = 8;
  nproc[2] = 1;

  res_x_total[0] = 128;
  res_x_total[1] = 128;
  res_x_total[2] = 1;

  for (int i = 0; i < nspecies; ++i) {
    species[i].res_v[0] = 32;
    species[i].res_v[1] = 32;
    species[i].res_v[2] = 32;
  }

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 2;
  dimensionality_v = 3;

  xb[0] = 0.;
  xb[1] = 0.;
  xb[2] = 0.;

  xe[0] = setup_var["L"];
  xe[1] = setup_var["L"];
  xe[2] = 1.;

  // A1, A2
  species[kElectron].vb[0] = -7.; // for m=0.04
  species[kElectron].vb[1] = -7.;
  species[kElectron].vb[2] = -7.;
  species[kElectron].ve[0] =  7.;
  species[kElectron].ve[1] =  7.;
  species[kElectron].ve[2] =  7.;
  species[kIon].vb[0] = -1.5;
  species[kIon].vb[1] = -1.5;
  species[kIon].vb[2] = -1.5;
  species[kIon].ve[0] =  1.5;
  species[kIon].ve[1] =  1.5;
  species[kIon].ve[2] =  1.5;

  // B1, B2
  //species[kElectron].vb[0] = -15.; // for m=0.04
  //species[kElectron].vb[1] = -15.;
  //species[kElectron].vb[2] = -15.;
  //species[kElectron].ve[0] =  15.;
  //species[kElectron].ve[1] =  15.;
  //species[kElectron].ve[2] =  15.;
  //species[kIon].vb[0] = -3.;
  //species[kIon].vb[1] = -3.;
  //species[kIon].vb[2] = -3.;
  //species[kIon].ve[0] =  3.;
  //species[kIon].ve[1] =  3.;
  //species[kIon].ve[2] =  3.;

  // time
  t_end = setup_var["L"] / setup_var["delta_u"];

  // numerics
  cfl = 0.1;
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


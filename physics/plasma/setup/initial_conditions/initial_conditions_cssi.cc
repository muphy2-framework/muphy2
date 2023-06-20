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

  // CSSI

  default_scheme = kF5eF5iM;

  output_directory = "./output/";
  noutputs_vtk = 150;
  noutputs_csv = 150;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 0.04; // Fujimoto & Sydora: 0.01
  species[kIon].m = 1.;
  species[kElectron].T0 = 1. / 12.;
  species[kIon].T0 = 5. / 12.;
  c0 = 20.; // Fujimoto & Sydora: 27.
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "cssi";
  setup_var["noise_level"] = 1.e-4;
  setup_var["T_bg"] = 1.e-3;

  // Fujimoto & Sydora setup (m_e = 0.01)
  //setup_var["n_bg"] = 1.;
  //setup_var["delta_e"] = .1;
  //setup_var["delta_i"] = .36;
  //setup_var["B0_e"] = 1.;
  //setup_var["B0_i"] = .76;
  //setup_var["V0_e"] = -10.;
  //setup_var["V0_i"] = 2.1;

  // adapted (m_e = 0.04)
  setup_var["n_bg"] = 0.25;
  setup_var["delta_e"] = .3;
  setup_var["delta_i"] = .2;
  setup_var["B0_e"] = .3;
  setup_var["B0_i"] = .1;
  setup_var["V0_e"] = -4.;
  setup_var["V0_i"] = 2.;

  // domain
  nproc[0] = 2;
  nproc[1] = 4;
  nproc[2] = 1;

  res_x_total[0] = 128;
  res_x_total[1] = 256;
  res_x_total[2] = 1;

  for (int i = 0; i < nspecies; ++i) {
    species[i].res_v[0] = 96;
    species[i].res_v[1] = 96;
    species[i].res_v[2] = 1;
  }

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 2;
  dimensionality_v = 2;

  xb[0] = 0.;
  xb[1] = -M_PI;
  xb[2] = 0.;

  xe[0] = M_PI;
  xe[1] = M_PI;
  xe[2] = 1.;

  // vspace should be larger and resolution higher
  species[kElectron].vb[0] = -12.;
  species[kElectron].vb[1] = -12.;
  species[kElectron].vb[2] = -12.;

  species[kElectron].ve[0] = 12.;
  species[kElectron].ve[1] = 12.;
  species[kElectron].ve[2] = 12.;

  species[kIon].vb[0] = -3.;
  species[kIon].vb[1] = -3.;
  species[kIon].vb[2] = -3.;

  species[kIon].ve[0] = 3.;
  species[kIon].ve[1] = 3.;
  species[kIon].ve[2] = 3.;

  // time
  t_end = 30.;

  // numerics
  cfl = 0.1;
  ten_moment_closure = "gradient_P";
  species[kElectron].ten_moment_k0 = 1./3. / sqrt(species[kElectron].m);
  species[kIon].ten_moment_k0 = 1./3. / sqrt(species[kIon].m);
  gradient_closure_subcycles = 8; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 15;

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


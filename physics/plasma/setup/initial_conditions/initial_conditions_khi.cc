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

  // Kelvin-Helmholtz Instability

  default_scheme = kF10eF10iM;

  output_directory = "./output/";
  noutputs_vtk = 150;
  noutputs_csv = 150;

  // physical parameters
  nspecies = 2;
  species = new Species[nspecies];
  species[kElectron].q = -1.;
  species[kIon].q = 1.;
  species[kElectron].m = 1.;
  species[kIon].m = 25.;
  c0 = 1.;
  mu0 = 1.;
  eps0 = 1./(c0*c0);

  // setup specific
  setup = "khi";
  setup_var["lambda"] = 10. * sqrt(species[kIon].m);
  setup_var["T0e_T0i"] = 0.2;
  setup_var["wpe_wce"] = 4.;

  // domain
  nproc[0] = 2;
  nproc[1] = 4;
  nproc[2] = 1;

  res_x_total[0] = 1600;
  res_x_total[1] = 2400;
  res_x_total[2] = 1;

  species[kElectron].res_v[0] = 32;
  species[kElectron].res_v[1] = 32;
  species[kElectron].res_v[2] = 1;
  species[kIon].res_v[0] = 24;
  species[kIon].res_v[1] = 24;
  species[kIon].res_v[2] = 1;

  bd[0] = 2;
  bd[1] = 2;
  bd[2] = 2;

  dimensionality_x = 2;
  dimensionality_v = 2;

  xb[0] = -9.*setup_var["lambda"];
  xb[1] = 0.;
  xb[2] = 0.;

  xe[0] = 9.*setup_var["lambda"];
  xe[1] = 27.*setup_var["lambda"];
  xe[2] = 1.;

  species[kElectron].vb[0] = -2.3;
  species[kElectron].vb[1] = -2.3;
  species[kElectron].vb[2] = -2.3;

  species[kElectron].ve[0] = 2.3;
  species[kElectron].ve[1] = 2.3;
  species[kElectron].ve[2] = 2.3;

  species[kIon].vb[0] = -1.;
  species[kIon].vb[1] = -1.;
  species[kIon].vb[2] = -1.;

  species[kIon].ve[0] = 1.;
  species[kIon].ve[1] = 1.;
  species[kIon].ve[2] = 1.;

  // time
  t_end = 300.*setup_var["wpe_wce"]*species[kIon].m; // 300 * \omega_p,e / \Omega_e * m_i/m_e

  // numerics
  cfl = 0.25;
  ten_moment_closure = "gradient_P";
  species[kElectron].ten_moment_k0 = 1./3. / sqrt(species[kElectron].m);
  species[kIon].ten_moment_k0 = 1./3. / sqrt(species[kIon].m);
  gradient_closure_subcycles = 8; // rule of thumb: double subcycles when doubling the resolution
  initial_maxwell_steps_per_fluid_step = 15;
  cweno_calc_j_from_flux = false; // makes boundary conditions easier

  // boundary conditions
  is_periodic[0] = false;
  is_periodic[1] = true;
  is_periodic[2] = true;

  bd_cond_v[0] = "zzpppp";
  bd_cond_v[1] = "sspppp";
  bd_cond_v[2] = "sspppp";

  bd_cond_E[0] = "sspppp";
  bd_cond_E[1] = "zzpppp";
  bd_cond_E[2] = "zzpppp";

  bd_cond_B[0] = "zzpppp";
  bd_cond_B[1] = "sspppp";
  bd_cond_B[2] = "sspppp";

  bd_cond_n[0] = "sspppp";
  bd_cond_nuu[0] = "sspppp";

  bd_cond_P[0] = "sspppp"; // xx
  bd_cond_P[1] = "zzpppp"; // xy
  bd_cond_P[2] = "zzpppp"; // xz
  bd_cond_P[3] = "sspppp"; // yy
  bd_cond_P[4] = "sspppp"; // yz
  bd_cond_P[5] = "sspppp"; // zz

  bd_cond_Q_h[0] = "zzpppp"; // xxx
  bd_cond_Q_h[1] = "sspppp"; // xxy
  bd_cond_Q_h[2] = "sspppp"; // xxz
  bd_cond_Q_h[3] = "zzpppp"; // xyy
  bd_cond_Q_h[4] = "zzpppp"; // xyz
  bd_cond_Q_h[5] = "zzpppp"; // xzz
  bd_cond_Q_h[6] = "sspppp"; // yyy
  bd_cond_Q_h[7] = "sspppp"; // yyz
  bd_cond_Q_h[8] = "sspppp"; // yzz
  bd_cond_Q_h[9] = "sspppp"; // zzz
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


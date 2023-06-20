/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "parameter.h"

Parameter::Parameter() {

  init();

  // domain related
  for (int i = 0; i < 3; ++i) {
    res_x[i] = res_x_total[i] / nproc[i];
    res_x_minus_one[i] = res_x[i] - 1;
    dx[i] = (xe[i] - xb[i]) / res_x_total[i];
    for (int s = 0; s < nspecies; ++s) {
      species[s].dv[i] = (species[s].ve[i] - species[s].vb[i]) / species[s].res_v[i];
    }
  }

  ncells_x_nobd = res_x[0] * res_x[1] * res_x[2];
  ncells_x = (res_x[0] + 2 * bd[0]) * (res_x[1] + 2 * bd[1]) * (res_x[2] + 2 * bd[2]);
  for (int s = 0; s < nspecies; ++s) {
    species[s].ncells_xv_nobd = ncells_x_nobd * species[s].res_v[0] * species[s].res_v[1] * species[s].res_v[2];
    species[s].ncells_xv = ncells_x * species[s].res_v[0] * species[s].res_v[1] * species[s].res_v[2];
  }

  // boundary conditions related
  // five moments
  bd_cond_five_moments[0] = bd_cond_n[0];
  for (int i = 0; i < 3; ++i)
    bd_cond_five_moments[1 + i] = bd_cond_v[i];
  bd_cond_five_moments[4] = bd_cond_nuu[0];

  // ten moments
  bd_cond_ten_moments[0] = bd_cond_n[0];
  for (int i = 0; i < 3; ++i)
    bd_cond_ten_moments[1 + i] = bd_cond_v[i];
  for (int i = 0; i < 6; ++i)
    bd_cond_ten_moments[4 + i] = bd_cond_P[i];

  // numerics related
  if (cweno_epsilon_deltaf < 0.) {
    cweno_epsilon_deltaf = cweno_epsilon;
  }
}

Parameter::~Parameter() {
  delete[] species;
  if (criterion != nullptr) {
    delete criterion;
  }
}


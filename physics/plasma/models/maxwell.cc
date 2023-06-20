/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "maxwell.h"
#include <math.h>
#include "framework/parameter.h"

namespace model {

Maxwell::Maxwell(const Parameter &parameter) : parameter_(parameter) {
  E = new real[3 * parameter_.ncells_x]();
  B = new real[3 * parameter_.ncells_x]();
  B_ipol = new real[3 * parameter_.ncells_x]();

  #pragma acc enter data copyin(E[0:3*parameter_.ncells_x])
  #pragma acc enter data copyin(B[0:3*parameter_.ncells_x])
  #pragma acc enter data copyin(B_ipol[0:3*parameter_.ncells_x])
}

Maxwell::~Maxwell() {
  #pragma acc exit data delete(E[0:3*parameter_.ncells_x])
  #pragma acc exit data delete(B[0:3*parameter_.ncells_x])
  #pragma acc exit data delete(B_ipol[0:3*parameter_.ncells_x])

  delete[] B_ipol;
  delete[] B;
  delete[] E;
}

int Maxwell::get_model_id() { return parameter_.kMaxwell; }

real Maxwell::get_max_dt() {
  const Parameter &p(parameter_);
  real sum_dx_sq = 0.;

  for (int i = 0; i < p.dimensionality_x; ++i)
    sum_dx_sq += p.dx[i] * p.dx[i];

  return p.cfl_maxwell / (p.c0 / sqrt(sum_dx_sq));
}

void Maxwell::step_E(const real dt, const real *j) {
  const Parameter &p(parameter_);
  step_e_maxwell_fdtd_(E, B, j, p.c0, p.mu0, p.res_x_minus_one, p.bd, p.dx, dt);
}

void Maxwell::step_B(const real dt) {
  const Parameter &p(parameter_);
  step_b_maxwell_fdtd_(E, B, p.res_x_minus_one, p.bd, p.dx, dt);
}

void Maxwell::step_B_ipol(const real dt) {
  const Parameter &p(parameter_);
  step_b_ipol_maxwell_fdtd_(E, B, B_ipol, p.res_x_minus_one, p.bd, p.dx, dt);
}

} // namespace model

/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "mhd.h"
#include "framework/parameter.h"

// TODO this model is work in progress and not functional

namespace model {

MHD1Temperature::MHD1Temperature(const Parameter &parameter)
  : parameter_(parameter) {
  five_moments = new real[5*parameter_.ncells_x]();
}

MHD1Temperature::~MHD1Temperature() { delete[] five_moments; }

int MHD1Temperature::get_model_id() { return parameter_.kMHD1Temperature; }

real MHD1Temperature::get_max_dt() {
  const Parameter &p(parameter_);
  return five_moment_mhd_max_dt_(five_moments, p.cfl, p.cweno_limiter,
                                 p.res_x_minus_one, p.bd, p.dx);
}

void MHD1Temperature::single_step_rk3(const real dt, real *five_moments_tmp,
                                      const real *E, const real *B, const real *J,
                                      const real alpha0, const real alpha1,
                                      const real beta) {
  const Parameter &p(parameter_);

  five_moment_mhd_single_step_rk3_(five_moments, five_moments_tmp, E, B, J, 
                                   alpha0, alpha1, beta, p.cweno_limiter,
                                   p.cweno_epsilon, p.res_x_minus_one, p.bd, p.dx,
                                   dt, p.dimensionality_x);
}

} // namespace model

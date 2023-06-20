/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "fluid5.h"
#include "framework/parameter.h"

namespace model {

Fluid5::Fluid5(const Parameter &parameter, const Species &species)
    : parameter_(parameter), species_(species) {
  five_moments = new real[5*parameter_.ncells_x]();
  #pragma acc enter data copyin(five_moments[0:5*parameter_.ncells_x])
}

Fluid5::~Fluid5() {
  #pragma acc exit data delete(five_moments[0:5*parameter_.ncells_x])
  delete[] five_moments;
}

int Fluid5::get_model_id() { return parameter_.kFluid5; }


real Fluid5::get_max_dt() {
  const Parameter &p(parameter_);

  real cfl = p.cfl;
  if (p.cfl_fluid_ions > 0. && species_.q > 0.) {
    cfl = p.cfl_fluid_ions;
  }

  return five_moment_max_dt_(five_moments, cfl, p.cweno_limiter,
                             p.res_x_minus_one, p.bd, p.dx, p.dimensionality_v);
}

void Fluid5::single_step_rk3(const real dt, real *five_moments_tmp,
                             const real *E, const real *B, const real alpha0,
                             const real alpha1, const real beta) {
  const Parameter &p(parameter_);
  const Species &s(species_);
  real q_m = s.q / s.m;

  five_moment_single_step_rk3_(five_moments, five_moments_tmp, E, B, q_m,
                               alpha0, alpha1, beta, p.cweno_limiter,
                               p.cweno_epsilon, p.res_x_minus_one, p.bd, p.dx,
                               dt, p.dimensionality_x, p.dimensionality_v);
}

void Fluid5::calc_un_from_flux(real *un) {
  const Parameter &p(parameter_);
  five_moment_calc_un_from_flux_(un, five_moments, p.cweno_limiter,
                                 p.cweno_epsilon, p.res_x_minus_one,
                                 p.bd, p.dimensionality_v);
}

} // namespace model

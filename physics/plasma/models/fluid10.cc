/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "fluid10.h"
#include "framework/parameter.h"

namespace model {

Fluid10::Fluid10(const Parameter &parameter, const Species &species)
    : parameter_(parameter), species_(species) {

  ten_moments = new real[10*parameter_.ncells_x]();
  #pragma acc enter data copyin(ten_moments[0:10*parameter_.ncells_x])
}


Fluid10::~Fluid10() {
  #pragma acc exit data delete(ten_moments[0:10*parameter_.ncells_x])
  delete[] ten_moments;
}


int Fluid10::get_model_id() { return parameter_.kFluid10; }


real Fluid10::get_max_dt() {
  const Parameter &p(parameter_);

  real cfl = p.cfl;
  if (p.cfl_fluid_ions > 0. && species_.q > 0.) {
    cfl = p.cfl_fluid_ions;
  }

  return ten_moment_max_dt_(ten_moments, cfl, p.cweno_limiter,
                            p.res_x_minus_one, p.bd, p.dx);
}


void Fluid10::calc_single_step_rk3(real* ten_moments_tmp, real* total_flux,
                                   real* source, const real* E, const real* B) {
  const Parameter &p(parameter_);

  real q_m = species_.q / species_.m;

  ten_moment_calc_single_step_rk3_(
      ten_moments_tmp, total_flux, source, E, B, q_m, p.cweno_limiter,
      p.cweno_epsilon, p.dimensionality_x, p.res_x_minus_one, p.bd, p.dx);
}


void Fluid10::calc_single_step_rk3(real* ten_moments_tmp, real* total_flux, real* source,
                                   const real* Q_raw, const real* E, const real* B) {
  const Parameter &p(parameter_);

  real q_m = species_.q / species_.m;

  ten_moment_deltaf_calc_single_step_rk3_(
      ten_moments_tmp, total_flux, source, Q_raw, E, B, q_m, p.cweno_limiter,
      p.cweno_epsilon_deltaf, p.dimensionality_x, p.res_x_minus_one, p.bd, p.dx);
}


void Fluid10::apply_single_step_rk3(const real dt, real* ten_moments_tmp,
                                    real* total_flux, real* source,
                                    const real alpha0, const real alpha1,
                                    const real beta) {
  const Parameter &p(parameter_);

  ten_moment_apply_single_step_rk3_(
      ten_moments, ten_moments_tmp, total_flux, source, alpha0, alpha1, beta,
      p.cweno_limiter, p.dimensionality_v, p.res_x_minus_one, p.bd, dt);
}


void Fluid10::calc_un_from_flux(real* un, bool deltaf_scheme/*=false*/) {
  const Parameter &p(parameter_);

  real cweno_epsilon = p.cweno_epsilon;
  if (deltaf_scheme) {
    cweno_epsilon = p.cweno_epsilon_deltaf;
  }
  ten_moment_calc_un_from_flux_(un, ten_moments, p.cweno_limiter,
                                cweno_epsilon, p.res_x_minus_one, p.bd);
}


void Fluid10::add_divheatflux_to_source(const real* divQ_h, real* source) {
  const Parameter &p(parameter_);
  ten_moment_add_divheatflux_to_source_(divQ_h, source, species_.m, p.res_x_minus_one);
}


void Fluid10::isotropization_heat_flux_closure(const real* ten_moments_tmp, real* source) {
  const Parameter &p(parameter_);

  ten_moment_isotropization_heat_flux_closure_(ten_moments_tmp, source,
                                         species_.ten_moment_k0, p.cweno_limiter,
                                         p.dimensionality_v, p.res_x_minus_one, p.bd);
}


void Fluid10::gradient_heat_flux_closure_cycle(const int iteration, const real dt,
                                               const real* ten_moments_tmp,
                                               real* source, real* P_m_tmp) {
  const Parameter &p(parameter_);

  if (p.ten_moment_closure == "gradient" || p.ten_moment_closure == "gradient_P") {
    ten_moment_gradient_p_heat_flux_closure_cycle_(
      ten_moments_tmp, source, P_m_tmp, species_.ten_moment_k0, iteration,
      p.gradient_closure_subcycles, p.cweno_limiter, p.res_x_minus_one, p.bd,
      p.dx, p.dimensionality_x, p.dimensionality_v, dt);
  } else if (p.ten_moment_closure == "gradient_T") {
    ten_moment_gradient_t_heat_flux_closure_cycle_(
      ten_moments_tmp, source, P_m_tmp, species_.ten_moment_k0, iteration,
      p.gradient_closure_subcycles, p.cweno_limiter, p.res_x_minus_one, p.bd,
      p.dx, p.dimensionality_x, dt);
  } else if (p.ten_moment_closure == "gradient_sym") {
    ten_moment_gradient_sym_heat_flux_closure_cycle_(
      ten_moments_tmp, source, P_m_tmp, species_.ten_moment_k0, iteration,
      p.gradient_closure_subcycles, p.cweno_limiter, p.res_x_minus_one, p.bd,
      p.dx, dt);
  }
}

} // namespace model

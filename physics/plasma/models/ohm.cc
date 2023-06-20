/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "ohm.h"
#include <math.h>
#include "framework/parameter.h"

// TODO this model is work in progress and not functional

namespace model {

Ohm::Ohm(const Parameter &parameter) : parameter_(parameter) {
  E = new real[3 * parameter_.ncells_x]();
  B = new real[3 * parameter_.ncells_x]();
  J = new real[3 * parameter_.ncells_x]();
}

Ohm::~Ohm() {
  delete[] J;
  delete[] B;
  delete[] E;
}

inline int Ohm::get_model_id(){ 
  return parameter_.kOhm;
}
 
real Ohm::get_max_dt(){
  return 1.; // TODO;
}

void Ohm::step_B(const real dt) {
  const Parameter &p(parameter_);
  ohm_step_b_(E, B, p.res_x_minus_one, p.bd, p.dx, dt);
}

void Ohm::calc_J() {
  const Parameter &p(parameter_);
  ohm_calc_j_(J, B, p.res_x_minus_one, p.bd, p.dx);
}

////////// E FIELD CALCULATION //////////

void OhmResistive::calc_E_fc(const real *u) {
  const Parameter &p(parameter_);
  ohm_calc_e_fc_ideal_(E, B, u, p.res_x_minus_one, p.bd);
  if (p.mhd_resistivity != 0.) ohm_calc_e_fc_resistivityterm_(E, J, p.mhd_resistivity, p.res_x_minus_one, p.bd);
}

void OhmHall::calc_E_fc(const real *u, const real *rho) {
  const Parameter &p(parameter_);
  ohm_calc_e_fc_ideal_(E, B, u, p.res_x_minus_one, p.bd);
  ohm_calc_e_fc_hallterm_(E, B, J, rho, p.res_x_minus_one, p.bd);
  if (p.mhd_resistivity != 0.) ohm_calc_e_fc_resistivityterm_(E, J, p.mhd_resistivity, p.res_x_minus_one, p.bd);
}

void OhmPressure::calc_E_fc(const real *u, const real *P, const real *rho) {
  const Parameter &p(parameter_);
  ohm_calc_e_fc_ideal_(E, B, u, p.res_x_minus_one, p.bd);
  ohm_calc_e_fc_pressureterm_(E, P, rho, p.dx, p.res_x_minus_one, p.bd);
  if (p.mhd_resistivity != 0.) ohm_calc_e_fc_resistivityterm_(E, J, p.mhd_resistivity, p.res_x_minus_one, p.bd);
}

void OhmHallPressure::calc_E_fc(const real *u, const real *P, const real *rho) {
  const Parameter &p(parameter_);
  ohm_calc_e_fc_ideal_(E, B, u, p.res_x_minus_one, p.bd);
  ohm_calc_e_fc_hallterm_(E, B, J, rho, p.res_x_minus_one, p.bd);
  ohm_calc_e_fc_pressureterm_(E, P, rho, p.dx, p.res_x_minus_one, p.bd);
  if (p.mhd_resistivity != 0.) ohm_calc_e_fc_resistivityterm_(E, J, p.mhd_resistivity, p.res_x_minus_one, p.bd);
}

} // namespace model

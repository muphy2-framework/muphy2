/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "MHD1TemperatureOhm.h"

#include <iostream>
#include <string.h>
#include <math.h>

#include "framework/interpolate.h"
#include "framework/mpi_boundary.h"
#include "framework/parameter.h"
#include "framework/block.h"
#include "physics/plasma/converters/converters.h"
#include "physics/plasma/models/mhd.h"
#include "physics/plasma/models/ohm.h"
#include "physics/plasma/setup/setup.h"
#include "physics/plasma/output/output.h"
#include "physics/plasma/criteria/select_scheme.h"

// TODO this scheme is work in progress and not functional

namespace scheme {
// one-temperature MHD, Ohm

MHD1TemperatureOhm::MHD1TemperatureOhm(const Block_id &block_id,
      Mpi_boundary &mpi_boundary, const Parameter &parameter)
    : block_id_(block_id), boundary_(mpi_boundary), parameter_(parameter) {

  mhd_ = new model::MHD1Temperature(parameter_);
  //ohm_ = new model::OhmResistive(parameter_);
  //ohm_ = new model::OhmHall(parameter_);
  ohm_ = new model::OhmPressure(parameter_);
  //ohm_ = new model::OhmHallPressure(parameter_);

  E_cc_ = new real[3 * parameter_.ncells_x]();
  B_cc_ = new real[3 * parameter_.ncells_x]();
  J_cc_ = new real[3 * parameter_.ncells_x]();

  five_moments_tmp_ = new real[5 * parameter_.ncells_x];
}


MHD1TemperatureOhm::~MHD1TemperatureOhm() {
  delete[] five_moments_tmp_;

  delete[] J_cc_;
  delete[] B_cc_;
  delete[] E_cc_;

  delete ohm_;
  delete mhd_;
}


int MHD1TemperatureOhm::get_scheme_id() { return parameter_.kMHD1TemperatureOhm; }


real MHD1TemperatureOhm::get_dt(int species) {
  return dt_;
}


void MHD1TemperatureOhm::init() {

  // load setup
  model::Model *plasma_models[1] = {mhd_};
  setup::init(plasma_models, ohm_, parameter_);

  // calc dt
  dt_ = std::min(mhd_->get_max_dt(), ohm_->get_max_dt());
  dt_ = boundary_.mpi_min(dt_);

  // calc j
  ohm_->calc_J();
  std::string orientation_J[3] = {"fcc", "cfc", "ccf"};
  boundary_.exchange_field(ohm_->J, parameter_.bd_cond_v, orientation_J, 3);
}


void MHD1TemperatureOhm::step() {

  // fluid
  std::string orientation_five_moments[5] = {"ccc", "ccc", "ccc", "ccc", "ccc"};

  ipol::yee_face_to_centered(ohm_->E, E_cc_, parameter_);
  ipol::yee_edge_to_centered(ohm_->B, B_cc_, parameter_);
  ipol::yee_face_to_centered(ohm_->J, J_cc_, parameter_);

  for (int i = 0; i < 5 * parameter_.ncells_x; ++i)
    five_moments_tmp_[i] = mhd_->five_moments[i];

  mhd_->single_step_rk3(dt_, five_moments_tmp_, E_cc_, B_cc_, J_cc_, 1., 0., 1.);
  boundary_.exchange_field(five_moments_tmp_, parameter_.bd_cond_five_moments, orientation_five_moments, 5);

  mhd_->single_step_rk3(dt_, five_moments_tmp_, E_cc_, B_cc_, J_cc_, .75, .25, .25);
  boundary_.exchange_field(five_moments_tmp_, parameter_.bd_cond_five_moments, orientation_five_moments, 5);

  mhd_->single_step_rk3(dt_, five_moments_tmp_, E_cc_, B_cc_, J_cc_, 1. / 3., 2. / 3., 2. / 3.);
  boundary_.exchange_field(five_moments_tmp_, parameter_.bd_cond_five_moments, orientation_five_moments, 5);

  for (int i = 0; i < 5 * parameter_.ncells_x; ++i)
    mhd_->five_moments[i] = five_moments_tmp_[i];

  // update dt
  if (!parameter_.freeze_dt) {
    dt_ = std::min(mhd_->get_max_dt(), ohm_->get_max_dt());
    dt_ = boundary_.mpi_min(dt_);
  }

  // ohm
  std::string orientation_E[3] = {"fcc", "cfc", "ccf"};
  std::string orientation_J[3] = {"fcc", "cfc", "ccf"};
  std::string orientation_B[3] = {"cff", "fcf", "ffc"};

  ohm_->step_B(dt_);  // dB/dt = - curl E
  boundary_.exchange_field(ohm_->B, parameter_.bd_cond_B, orientation_B, 3);

  ohm_->calc_J();  // J = curl B
  boundary_.exchange_field(ohm_->J, parameter_.bd_cond_v, orientation_J, 3);

  real *u = new real[3*parameter_.ncells_x];
  real *P_e = new real[3 * parameter_.ncells_x];
  convert::un_to_u(u, &mhd_->five_moments[parameter_.ncells_x], &mhd_->five_moments[0], parameter_);
  convert::mhd_1temperature_to_electron_pressure(P_e, mhd_->five_moments, parameter_);
  //ohm_->calc_E_fc(u); // E = - u x B, with ipol of u,B to fc
  //ohm_->calc_E_fc(u, &mhd_->five_moments[0]); // E = - u x B, with ipol of u,B to fc
  ohm_->calc_E_fc(u, P_e, &mhd_->five_moments[0]); // E = - u x B, with ipol of u,B to fc
  boundary_.exchange_field(ohm_->E, parameter_.bd_cond_E, orientation_E, 3);

  delete[] P_e;
  delete[] u;
}


void MHD1TemperatureOhm::output(real t, int output_number, bool output_vtk) {
  model::Model *plasma_models[2] = {mhd_, mhd_};
  output::prepare_and_write(plasma_models, ohm_, 2, t, output_number, output_vtk, block_id_, parameter_);
}


int MHD1TemperatureOhm::evaluate_criterion() {
  model::Model *plasma_models[2] = {mhd_, mhd_};
  return criterion::evaluate_criterion(plasma_models, ohm_, block_id_, parameter_);
}

} // namespace scheme

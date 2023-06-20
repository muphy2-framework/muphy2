/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "VeViPdeltaf.h"

#include <iostream>
#include <math.h>
#include <string.h>

#include "framework/utilities.h"
#include "framework/restart.h"
#include "framework/interpolate.h"
#include "framework/mpi_boundary.h"
#include "framework/parameter.h"
#include "framework/block.h"
#include "physics/plasma/converters/converters.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/models/fluid5.h"
#include "physics/plasma/models/poisson.h"
#include "physics/plasma/setup/setup.h"
#include "physics/plasma/output/output.h"
#include "physics/plasma/criteria/select_scheme.h"

namespace scheme {

// Vlasov electrons, Vlasov ions, Poisson

VeViPdeltaf::VeViPdeltaf(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter)
    : block_id_(block_id), boundary_(mpi_boundary), parameter_(parameter) {

  vlasov_e_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kElectron]);
  vlasov_i_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kIon]);
  fluid_e_  = new model::Fluid10(parameter_, parameter_.species[parameter_.kElectron]);
  fluid_i_  = new model::Fluid10(parameter_, parameter_.species[parameter_.kIon]);
  poisson_ = new model::Poisson(parameter_);

  rho_ = new real[parameter_.ncells_x]();
  E_old_ = new real[3 * parameter_.ncells_x]();

  Q_raw_old_e_ = new real[10 * parameter_.ncells_x]();
  Q_raw_old_i_ = new real[10 * parameter_.ncells_x]();
  Q_raw_new_e_ = new real[10 * parameter_.ncells_x]();
  Q_raw_new_i_ = new real[10 * parameter_.ncells_x]();

  ten_moments_tmp_e_ = new real[10 * parameter_.ncells_x]();
  ten_moments_tmp_i_ = new real[10 * parameter_.ncells_x]();
  total_flux_e_ = new real[10 * parameter_.ncells_x_nobd]();
  total_flux_i_ = new real[10 * parameter_.ncells_x_nobd]();
  source_e_ = new real[10 * parameter_.ncells_x_nobd]();
  source_i_ = new real[10 * parameter_.ncells_x_nobd]();

  #pragma acc enter data copyin(rho_[0:parameter_.ncells_x], \
                                E_old_[0:3*parameter_.ncells_x], \
                                Q_raw_old_e_[0:10*parameter_.ncells_x], \
                                Q_raw_old_i_[0:10*parameter_.ncells_x], \
                                Q_raw_new_e_[0:10*parameter_.ncells_x], \
                                Q_raw_new_i_[0:10*parameter_.ncells_x], \
                                ten_moments_tmp_e_[0:10*parameter_.ncells_x], \
                                ten_moments_tmp_i_[0:10*parameter_.ncells_x], \
                                total_flux_e_[0:10*parameter_.ncells_x_nobd], \
                                total_flux_i_[0:10*parameter_.ncells_x_nobd], \
                                source_e_[0:10*parameter_.ncells_x_nobd], \
                                source_i_[0:10*parameter_.ncells_x_nobd])
}


VeViPdeltaf::~VeViPdeltaf() {
  #pragma acc exit data delete(rho_[0:parameter_.ncells_x], \
                               E_old_[0:3*parameter_.ncells_x], \
                               Q_raw_old_e_[0:10*parameter_.ncells_x], \
                               Q_raw_old_i_[0:10*parameter_.ncells_x], \
                               Q_raw_new_e_[0:10*parameter_.ncells_x], \
                               Q_raw_new_i_[0:10*parameter_.ncells_x], \
                               ten_moments_tmp_e_[0:10*parameter_.ncells_x], \
                               ten_moments_tmp_i_[0:10*parameter_.ncells_x], \
                               total_flux_e_[0:10*parameter_.ncells_x_nobd], \
                               total_flux_i_[0:10*parameter_.ncells_x_nobd], \
                               source_e_[0:10*parameter_.ncells_x_nobd], \
                               source_i_[0:10*parameter_.ncells_x_nobd])

  delete[] source_i_;
  delete[] source_e_;
  delete[] total_flux_i_;
  delete[] total_flux_e_;
  delete[] ten_moments_tmp_i_;
  delete[] ten_moments_tmp_e_;

  delete[] Q_raw_new_i_;
  delete[] Q_raw_new_e_;
  delete[] Q_raw_old_i_;
  delete[] Q_raw_old_e_;

  delete[] E_old_;
  delete[] rho_;

  delete poisson_;
  delete fluid_i_;
  delete fluid_e_;
  delete vlasov_i_;
  delete vlasov_e_;
}


int VeViPdeltaf::get_scheme_id() { return parameter_.kVeViPdeltaf; }


real VeViPdeltaf::get_dt(int species) { return dt_; }


void VeViPdeltaf::init() {

  if (parameter_.restart) {

    restart::load_field(vlasov_e_->f, 0, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(vlasov_i_->f, 1, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(fluid_e_->ten_moments, 2, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(fluid_i_->ten_moments, 3, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(poisson_->phi, 4, parameter_.output_directory, block_id_.my_coords);

    #pragma acc update device(vlasov_e_->f[0:parameter_.species[parameter_.kElectron].ncells_xv], \
                              vlasov_i_->f[0:parameter_.species[parameter_.kIon].ncells_xv], \
                              fluid_e_->ten_moments[0:10*parameter_.ncells_x], \
                              fluid_i_->ten_moments[0:10*parameter_.ncells_x], \
                              poisson_->phi[0:parameter_.ncells_x])
  } else {
    // load setup
    model::Model *plasma_models[2] =  {vlasov_e_, vlasov_i_};
    setup::init(plasma_models, poisson_, parameter_);
    convert::f_to_ten_moments(fluid_e_->ten_moments, vlasov_e_->f, parameter_.species[parameter_.kElectron], parameter_);
    convert::f_to_ten_moments(fluid_i_->ten_moments, vlasov_i_->f, parameter_.species[parameter_.kIon], parameter_);

    // calc rho
    convert::f_to_rho(rho_, vlasov_e_->f, vlasov_i_->f, parameter_);
    std::string orientation_rho[1] = {"ccc"};
    boundary_.exchange_field(rho_, parameter_.bd_cond_n, orientation_rho, 1);

    // poisson
    for (int i = 0; i < parameter_.poisson_schwarz_iterations; ++i) {

      poisson_->solver_iterations(rho_);
      real phi_average = poisson_->sum_phi() / parameter_.ncells_x_nobd;
      phi_average = boundary_.mpi_sum(phi_average) /
                    (parameter_.nproc[0]*parameter_.nproc[1]*parameter_.nproc[2]);
      util::array_add_scalar(poisson_->phi, -phi_average, parameter_.ncells_x);
      std::string orientation_phi[1] = {"ccc"};
      boundary_.exchange_schwarz(poisson_->phi, parameter_.bd_cond_phi, orientation_phi);
    }
  }

  poisson_->calc_E_from_phi();

  // needed for Vlasov's first half-step v
  std::string orientation_E[3] = {"ccc", "ccc", "ccc"};
  boundary_.exchange_field(poisson_->E, parameter_.bd_cond_E, orientation_E, 3);

  util::copy_array(ten_moments_tmp_e_, fluid_e_->ten_moments, 10*parameter_.ncells_x);
  util::copy_array(ten_moments_tmp_i_, fluid_i_->ten_moments, 10*parameter_.ncells_x);

  // calc heat flux
  std::string orientation_Q_h[10] = {"ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc"};
  bool exchange_in_this_direction_e[6] = {0, 0, 0, 0, 0, 0};
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
      exchange_in_this_direction_e[d] = true;
    }
  }
  convert::f_to_raw_heatflux(Q_raw_new_e_, vlasov_e_->f, parameter_.species[parameter_.kElectron],
                             exchange_in_this_direction_e, parameter_);
  boundary_.exchange_field(Q_raw_new_e_, parameter_.bd_cond_Q_h, orientation_Q_h, exchange_in_this_direction_e, 10);

  bool exchange_in_this_direction_i[6] = {0, 0, 0, 0, 0, 0};
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
      exchange_in_this_direction_i[d] = true;
    }
  }
  convert::f_to_raw_heatflux(Q_raw_new_i_, vlasov_i_->f, parameter_.species[parameter_.kIon],
                             exchange_in_this_direction_i, parameter_);
  boundary_.exchange_field(Q_raw_new_i_, parameter_.bd_cond_Q_h, orientation_Q_h, exchange_in_this_direction_i, 10);


  // calc dt
  real dt_vlasov = std::min(vlasov_e_->get_max_dt(), vlasov_i_->get_max_dt());
  real dt_fluid = std::min(fluid_e_->get_max_dt(), fluid_i_->get_max_dt());
  dt_ = std::min(dt_vlasov, dt_fluid);
}


void VeViPdeltaf::step() {

  // for time interpolation
  util::copy_array(E_old_, poisson_->E, 3*parameter_.ncells_x);
  util::copy_array(Q_raw_old_e_, Q_raw_new_e_, 10*parameter_.ncells_x);
  util::copy_array(Q_raw_old_i_, Q_raw_new_i_, 10*parameter_.ncells_x);


  // vlasov half step v
  if (parameter_.vlasov_solver == "pfc") { // backsubstitution method needs no additional splitting
    vlasov_e_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    vlasov_i_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    if (parameter_.dimensionality_v > 1) {
      vlasov_e_->step_vy(dt_*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
      vlasov_i_->step_vy(dt_*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    }
    if (parameter_.dimensionality_v > 2) {
      vlasov_e_->step_vz(dt_*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
      vlasov_i_->step_vz(dt_*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    }
  }
  else { // splitting for other solvers
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
      vlasov_i_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    }
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
      vlasov_i_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    }
    vlasov_e_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B, include_boundary_cells_); // 1. x
    vlasov_i_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
      vlasov_i_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    }
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
      vlasov_i_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B, include_boundary_cells_);
    }
  }

  // vlasov step x
  vlasov_e_->step_x(dt_);
  vlasov_i_->step_x(dt_);
  if (parameter_.dimensionality_x > 1) {
    vlasov_e_->step_y(dt_);
    vlasov_i_->step_y(dt_);
  }
  if (parameter_.dimensionality_x > 2) {
    vlasov_e_->step_z(dt_);
    vlasov_i_->step_z(dt_);
  }


  // calc rho
  convert::f_to_rho(rho_, vlasov_e_->f, vlasov_i_->f, parameter_);
  std::string orientation_rho[1] = {"ccc"};
  boundary_.exchange_field(rho_, parameter_.bd_cond_n, orientation_rho, 1);

  // poisson
  for (int i = 0; i < parameter_.poisson_schwarz_iterations; ++i) {

    poisson_->solver_iterations(rho_);
    real phi_average = poisson_->sum_phi() / parameter_.ncells_x_nobd;
    phi_average = boundary_.mpi_sum(phi_average) /
                  (parameter_.nproc[0]*parameter_.nproc[1]*parameter_.nproc[2]);
    util::array_add_scalar(poisson_->phi, -phi_average, parameter_.ncells_x);
    std::string orientation_phi[1] = {"ccc"};
    boundary_.exchange_schwarz(poisson_->phi, parameter_.bd_cond_phi, orientation_phi);
  }
  poisson_->calc_E_from_phi();

  // needed for Vlasov's first half-step v
  std::string orientation_E[3] = {"ccc", "ccc", "ccc"};
  boundary_.exchange_field(poisson_->E, parameter_.bd_cond_E, orientation_E, 3);
 

  // vlasov second half step v
  if (parameter_.vlasov_solver == "pfc") { // backsubstitution method needs no additional splitting
    vlasov_e_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B);
    vlasov_i_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B);
    if (parameter_.dimensionality_v > 1) {
      vlasov_e_->step_vy(dt_*0.5, *poisson_->E, *poisson_->B);
      vlasov_i_->step_vy(dt_*0.5, *poisson_->E, *poisson_->B);
    }
    if (parameter_.dimensionality_v > 2) {
      vlasov_e_->step_vz(dt_*0.5, *poisson_->E, *poisson_->B);
      vlasov_i_->step_vz(dt_*0.5, *poisson_->E, *poisson_->B);
    }
  }
  else { // splitting for other solvers
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
      vlasov_i_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
    }
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
      vlasov_i_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
    }
    vlasov_e_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B); // 1. x
    vlasov_i_->step_vx(dt_ * 0.5, *poisson_->E, *poisson_->B);
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
      vlasov_i_->step_vy(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
    }
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
      vlasov_i_->step_vz(dt_*0.5*0.5, *poisson_->E, *poisson_->B);
    }
  }


  // calc heat flux
  std::string orientation_Q_h[10] = {"ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc"};
  bool exchange_in_this_direction_e[6] = {0, 0, 0, 0, 0, 0};
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
      exchange_in_this_direction_e[d] = true;
    }
  }
  convert::f_to_raw_heatflux(Q_raw_new_e_, vlasov_e_->f, parameter_.species[parameter_.kElectron],
                             exchange_in_this_direction_e, parameter_);
  boundary_.exchange_field(Q_raw_new_e_, parameter_.bd_cond_Q_h, orientation_Q_h, exchange_in_this_direction_e, 10);

  bool exchange_in_this_direction_i[6] = {0, 0, 0, 0, 0, 0};
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
      exchange_in_this_direction_i[d] = true;
    }
  }
  convert::f_to_raw_heatflux(Q_raw_new_i_, vlasov_i_->f, parameter_.species[parameter_.kIon],
                             exchange_in_this_direction_i, parameter_);
  boundary_.exchange_field(Q_raw_new_i_, parameter_.bd_cond_Q_h, orientation_Q_h, exchange_in_this_direction_i, 10);


  // fluid rk step 1
  fluid_e_->calc_single_step_rk3(ten_moments_tmp_e_, total_flux_e_, source_e_, Q_raw_old_e_, E_old_, poisson_->B);
  fluid_e_->apply_single_step_rk3(dt_, ten_moments_tmp_e_, total_flux_e_, source_e_, 1., 0., 1.);
  fluid_i_->calc_single_step_rk3(ten_moments_tmp_i_, total_flux_i_, source_i_, Q_raw_old_i_, E_old_, poisson_->B);
  fluid_i_->apply_single_step_rk3(dt_, ten_moments_tmp_i_, total_flux_i_, source_i_, 1., 0., 1.);

  exchange_plasma_solvers(vlasov_e_->f, vlasov_i_->f, ten_moments_tmp_e_, ten_moments_tmp_i_, runge_kutta_only_);

  // fluid rk step 2
  fluid_e_->calc_single_step_rk3(ten_moments_tmp_e_, total_flux_e_, source_e_, Q_raw_new_e_, poisson_->E, poisson_->B);
  fluid_e_->apply_single_step_rk3(dt_, ten_moments_tmp_e_, total_flux_e_, source_e_, .75, .25, .25);
  fluid_i_->calc_single_step_rk3(ten_moments_tmp_i_, total_flux_i_, source_i_, Q_raw_new_i_, poisson_->E, poisson_->B);
  fluid_i_->apply_single_step_rk3(dt_, ten_moments_tmp_i_, total_flux_i_, source_i_, .75, .25, .25);

  exchange_plasma_solvers(vlasov_e_->f, vlasov_i_->f, ten_moments_tmp_e_, ten_moments_tmp_i_, runge_kutta_only_);

  // interpolate E and heat flux to intermediate time: E_old = 0.5*E_old + 0.5*E_new, Q_old = 0.5*Q_old + 0.5*Q_new
  util::add_arrays(E_old_, E_old_, poisson_->E, 0.5, 0.5, 3*parameter_.ncells_x);
  util::add_arrays(Q_raw_old_e_, Q_raw_old_e_, Q_raw_new_e_, 0.5, 0.5, 10*parameter_.ncells_x);
  util::add_arrays(Q_raw_old_i_, Q_raw_old_i_, Q_raw_new_i_, 0.5, 0.5, 10*parameter_.ncells_x);

  // fluid rk step 3
  fluid_e_->calc_single_step_rk3(ten_moments_tmp_e_, total_flux_e_, source_e_, Q_raw_old_e_, E_old_, poisson_->B);
  fluid_e_->apply_single_step_rk3(dt_, ten_moments_tmp_e_, total_flux_e_, source_e_, 1./3., 2./3., 2./3.);
  fluid_i_->calc_single_step_rk3(ten_moments_tmp_i_, total_flux_i_, source_i_, Q_raw_old_i_, E_old_, poisson_->B);
  fluid_i_->apply_single_step_rk3(dt_, ten_moments_tmp_i_, total_flux_i_, source_i_, 1./3., 2./3., 2./3.);


  convert::delta_f_correction(vlasov_e_->f, ten_moments_tmp_e_,
                              parameter_.species[parameter_.kElectron], parameter_);
  convert::delta_f_correction(vlasov_i_->f, ten_moments_tmp_i_,
                              parameter_.species[parameter_.kIon], parameter_);

  exchange_plasma_solvers(vlasov_e_->f, vlasov_i_->f, ten_moments_tmp_e_, ten_moments_tmp_i_);

  util::copy_array(fluid_e_->ten_moments, ten_moments_tmp_e_, 10*parameter_.ncells_x);
  util::copy_array(fluid_i_->ten_moments, ten_moments_tmp_i_, 10*parameter_.ncells_x);

  // update dt
  if (!parameter_.freeze_dt) {
    real dt_vlasov = std::min(vlasov_e_->get_max_dt(), vlasov_i_->get_max_dt());
    real dt_fluid = std::min(fluid_e_->get_max_dt(), fluid_i_->get_max_dt());
    dt_ = std::min(dt_vlasov, dt_fluid);
  }
}


void VeViPdeltaf::exchange_plasma_solvers(real* data_vlasov_e, real* data_vlasov_i,
        real* data_fluid_e, real* data_fluid_i, bool runge_kutta_only) {

  real* buffer_out_e[parameter_.kNDirections],
      * buffer_out_i[parameter_.kNDirections],
      * buffer_in_e [parameter_.kNDirections],
      * buffer_in_i [parameter_.kNDirections];
  int buffer_size_e[parameter_.kNDirections],
      buffer_size_i[parameter_.kNDirections];

  // fluid exchange
  std::string orientation[10] = {"ccc","ccc","ccc","ccc","ccc","ccc","ccc","ccc","ccc","ccc"};

  bool exchange_this_direction[6] = {0, 0, 0, 0, 0, 0};

  for (int d = 0; d < parameter_.kNDirections; ++d) {

    // prepare buffers
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {

      exchange_this_direction[d] = true;
      // electrons
      buffer_size_e[d] = boundary_.calc_buffer_size(d, 10);
      buffer_out_e[d] = new real[buffer_size_e[d]];
      buffer_in_e[d]  = new real[buffer_size_e[d]];
      #pragma acc enter data create(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                    buffer_in_e [d:1][0:buffer_size_e[d]])
      boundary_.extract_buffer(data_fluid_e, buffer_out_e[d], d, 10);
      // ions
      buffer_size_i[d] = boundary_.calc_buffer_size(d, 10);
      buffer_out_i[d] = new real[buffer_size_i[d]];
      buffer_in_i[d]  = new real[buffer_size_i[d]];
      #pragma acc enter data create(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                    buffer_in_i [d:1][0:buffer_size_i[d]])
      boundary_.extract_buffer(data_fluid_i, buffer_out_i[d], d, 10);
    }
  }

  // send/receive buffers
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (exchange_this_direction[d])
      boundary_.exchange_buffer(buffer_out_e[d], buffer_in_e[d], buffer_size_e[d], d);
  }
  boundary_.wait(exchange_this_direction);
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (exchange_this_direction[d])
      boundary_.exchange_buffer(buffer_out_i[d], buffer_in_i[d], buffer_size_i[d], d);
  }
  boundary_.wait(exchange_this_direction);

  for (int d = 0; d < parameter_.kNDirections; ++d) {

    // apply received buffers (convert if necessary)
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
      // electrons
      boundary_.apply_buffer(data_fluid_e, buffer_in_e[d], d, 10);
      // ions
      boundary_.apply_buffer(data_fluid_i, buffer_in_i[d], d, 10);
    }

    // clean up
    if (exchange_this_direction[d]) {
      #pragma acc exit data delete(buffer_out_e[d:1][0:buffer_size_e[d]], buffer_out_i[d:1][0:buffer_size_i[d]],\
                                   buffer_in_e [d:1][0:buffer_size_e[d]], buffer_in_i [d:1][0:buffer_size_i[d]])
      delete[] buffer_out_e[d];
      delete[] buffer_in_e [d];
      delete[] buffer_out_i[d];
      delete[] buffer_in_i [d];
    }

    // apply boundary conditions
    boundary_.apply_boundary_conditions(data_fluid_e, parameter_.bd_cond_ten_moments,
            orientation, d, 10);
    boundary_.apply_boundary_conditions(data_fluid_i, parameter_.bd_cond_ten_moments,
            orientation, d, 10);
  }

  // Vlasov exchange
  if (!runge_kutta_only) {
    std::string orientation[1] = {"ccc"};

    bool exchange_this_direction_e[6] = {0, 0, 0, 0, 0, 0},
         exchange_this_direction_i[6] = {0, 0, 0, 0, 0, 0};

    for (int d = 0; d < parameter_.kNDirections; ++d) {

      // prepare buffers
      if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
        // electrons
        exchange_this_direction_e[d] = true;
        buffer_size_e[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kElectron);
        buffer_out_e[d] = new real[buffer_size_e[d]];
        buffer_in_e[d]  = new real[buffer_size_e[d]];
        #pragma acc enter data create(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                      buffer_in_e [d:1][0:buffer_size_e[d]])
        boundary_.extract_buffer(data_vlasov_e, buffer_out_e[d], d, 1, true, parameter_.kElectron);
        // ions
        exchange_this_direction_i[d] = true;
        buffer_size_i[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kIon);
        buffer_out_i[d] = new real[buffer_size_i[d]];
        buffer_in_i[d]  = new real[buffer_size_i[d]];
        #pragma acc enter data create(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                      buffer_in_i [d:1][0:buffer_size_i[d]])
        boundary_.extract_buffer(data_vlasov_i, buffer_out_i[d], d, 1, true, parameter_.kIon);
      }
    }

    // send/receive buffers
    for (int d = 0; d < parameter_.kNDirections; ++d) {
      if (exchange_this_direction_e[d])
        boundary_.exchange_buffer(buffer_out_e[d], buffer_in_e[d], buffer_size_e[d], d);
    }
    boundary_.wait(exchange_this_direction_e);
    for (int d = 0; d < parameter_.kNDirections; ++d) {
      if (exchange_this_direction_i[d])
        boundary_.exchange_buffer(buffer_out_i[d], buffer_in_i[d], buffer_size_i[d], d);
    }
    boundary_.wait(exchange_this_direction_i);

    for (int d = 0; d < parameter_.kNDirections; ++d) {

      if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViPdeltaf) {
        // electrons
        boundary_.apply_buffer(data_vlasov_e, buffer_in_e[d], d, 1, true, parameter_.kElectron);
        // ions
        boundary_.apply_buffer(data_vlasov_i, buffer_in_i[d], d, 1, true, parameter_.kIon);
      }

      // clean up
      if (exchange_this_direction_e[d]) {
        #pragma acc exit data delete(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                     buffer_in_e [d:1][0:buffer_size_e[d]])
        delete[] buffer_out_e[d];
        delete[] buffer_in_e [d];
      }
      if (exchange_this_direction_i[d]) {
        #pragma acc exit data delete(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                     buffer_in_i [d:1][0:buffer_size_i[d]])
        delete[] buffer_out_i[d];
        delete[] buffer_in_i [d];
      }

      // apply boundary conditions
      boundary_.apply_boundary_conditions(data_vlasov_e, parameter_.bd_cond_v,
              orientation, d, 1, true, parameter_.kElectron);
      boundary_.apply_boundary_conditions(data_vlasov_i, parameter_.bd_cond_v,
              orientation, d, 1, true, parameter_.kIon);
    }
  } // end if(!runge_kutta_only)
}


void VeViPdeltaf::output(real t, int output_number, bool output_vtk) {
  model::Model *plasma_models[4] = {vlasov_e_, vlasov_i_, fluid_e_, fluid_i_};
  output::prepare_and_write(plasma_models, poisson_, 4, t, output_number, output_vtk, block_id_, parameter_);
}


int VeViPdeltaf::evaluate_criterion() {
  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  return criterion::evaluate_criterion(plasma_models, poisson_, block_id_, parameter_);
}

} // namespace scheme

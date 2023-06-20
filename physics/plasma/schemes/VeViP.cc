/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "VeViP.h"

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
#include "physics/plasma/models/poisson.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/setup/setup.h"
#include "physics/plasma/output/output.h"
#include "physics/plasma/criteria/select_scheme.h"

namespace scheme {

// Vlasov electrons, Vlasov ions, Poisson

VeViP::VeViP(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter)
    : block_id_(block_id), boundary_(mpi_boundary), parameter_(parameter) {

  vlasov_e_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kElectron]);
  vlasov_i_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kIon]);
  poisson_ = new model::Poisson(parameter_);

  rho_ = new real[parameter_.ncells_x]();
  #pragma acc enter data copyin(rho_[0:parameter_.ncells_x])
}


VeViP::~VeViP() {
  #pragma acc exit data delete(rho_[0:parameter_.ncells_x])
  delete[] rho_;

  delete poisson_;
  delete vlasov_i_;
  delete vlasov_e_;
}


int VeViP::get_scheme_id() { return parameter_.kVeViP; }


real VeViP::get_dt(int species) { return dt_; }


void VeViP::init() {

  if (parameter_.restart) {

    restart::load_field(vlasov_e_->f, 0, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(vlasov_i_->f, 1, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(poisson_->phi, 2, parameter_.output_directory, block_id_.my_coords);
    #pragma acc update device(vlasov_e_->f[0:parameter_.species[parameter_.kElectron].ncells_xv], \
                              vlasov_i_->f[0:parameter_.species[parameter_.kIon].ncells_xv], \
                              poisson_->phi[0:parameter_.ncells_x])
  } else {
    // load setup
    model::Model *plasma_models[2] =  {vlasov_e_, vlasov_i_};
    setup::init(plasma_models, poisson_, parameter_);

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
  std::string orientation_E[3] = {"ccc","ccc","ccc"};
  boundary_.exchange_field(poisson_->E, parameter_.bd_cond_E, orientation_E, 3);

  // calc dt
  dt_ = std::min(vlasov_e_->get_max_dt(), vlasov_i_->get_max_dt());
  dt_ = boundary_.mpi_min(dt_);
}


void VeViP::step() {

  // vlasov
  // half step v
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


  // step x
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
  std::string orientation_E[3] = {"ccc","ccc","ccc"};
  boundary_.exchange_field(poisson_->E, parameter_.bd_cond_E, orientation_E, 3);


  // second half step v
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

  exchange_plasma_solvers(vlasov_e_->f, vlasov_i_->f);

  // update dt
  if (!parameter_.freeze_dt) {
    dt_ = std::min(vlasov_e_->get_max_dt(), vlasov_i_->get_max_dt());
    dt_ = boundary_.mpi_min(dt_);
  }
}


void VeViP::exchange_plasma_solvers(real* data_e, real* data_i) {

  real* buffer_out_e[parameter_.kNDirections],
      * buffer_out_i[parameter_.kNDirections],
      * buffer_in_e [parameter_.kNDirections],
      * buffer_in_i [parameter_.kNDirections];
  int buffer_size_e[parameter_.kNDirections],
      buffer_size_i[parameter_.kNDirections];

  std::string orientation[1] = {"ccc"};

  bool exchange_this_direction[6] = {0, 0, 0, 0, 0, 0};

  for (int d = 0; d < parameter_.kNDirections; ++d) {

    // prepare buffers
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViP) {
      exchange_this_direction[d] = true;
      // electrons
      buffer_size_e[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kElectron);
      buffer_out_e[d] = new real[buffer_size_e[d]];
      buffer_in_e[d]  = new real[buffer_size_e[d]];
      #pragma acc enter data create(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                    buffer_in_e [d:1][0:buffer_size_e[d]])
      boundary_.extract_buffer(data_e, buffer_out_e[d], d, 1, true, parameter_.kElectron);
      // ions
      buffer_size_i[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kIon);
      buffer_out_i[d] = new real[buffer_size_i[d]];
      buffer_in_i[d]  = new real[buffer_size_i[d]];
      #pragma acc enter data create(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                    buffer_in_i [d:1][0:buffer_size_i[d]])
      boundary_.extract_buffer(data_i, buffer_out_i[d], d, 1, true, parameter_.kIon);
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
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViP) {
      // electrons
      boundary_.apply_buffer(data_e, buffer_in_e[d], d, 1, true, parameter_.kElectron);
      // ions
      boundary_.apply_buffer(data_i, buffer_in_i[d], d, 1, true, parameter_.kIon);
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
    boundary_.apply_boundary_conditions(data_e, parameter_.bd_cond_v,
            orientation, d, 1, true, parameter_.kElectron);
    boundary_.apply_boundary_conditions(data_i, parameter_.bd_cond_v,
            orientation, d, 1, true, parameter_.kIon);
  }
}


void VeViP::output(real t, int output_number, bool output_vtk) {
  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  output::prepare_and_write(plasma_models, poisson_, 2, t, output_number, output_vtk, block_id_, parameter_);
}


int VeViP::evaluate_criterion() {
  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  return criterion::evaluate_criterion(plasma_models, poisson_, block_id_, parameter_);
}

} // namespace scheme

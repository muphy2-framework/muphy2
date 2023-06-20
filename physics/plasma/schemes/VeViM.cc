/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "VeViM.h"

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
#include "physics/plasma/models/maxwell.h"
#include "physics/plasma/setup/setup.h"
#include "physics/plasma/output/output.h"
#include "physics/plasma/criteria/select_scheme.h"
#include "physics/plasma/criteria/load_balancing.h"

namespace scheme {

// Vlasov electrons, Vlasov ions, Maxwell

VeViM::VeViM(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter)
    : block_id_(block_id), boundary_(mpi_boundary), parameter_(parameter) {

  vlasov_e_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kElectron]);
  vlasov_i_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kIon]);
  maxwell_  = new model::Maxwell(parameter_);

  j_        = new real[3 * parameter_.ncells_x]();
  j_fc_     = new real[3 * parameter_.ncells_x]();
  E_cc_     = new real[3 * parameter_.ncells_x]();
  B_cc_     = new real[3 * parameter_.ncells_x]();
  #pragma acc enter data copyin(j_   [0:3*parameter_.ncells_x], \
                                j_fc_[0:3*parameter_.ncells_x], \
                                E_cc_[0:3*parameter_.ncells_x], \
                                B_cc_[0:3*parameter_.ncells_x])
}


VeViM::~VeViM() {
  #pragma acc exit data delete(B_cc_[0:3*parameter_.ncells_x], \
                               E_cc_[0:3*parameter_.ncells_x], \
                               j_fc_[0:3*parameter_.ncells_x], \
                               j_   [0:3*parameter_.ncells_x])
  delete[] B_cc_;
  delete[] E_cc_;
  delete[] j_fc_;
  delete[] j_;

  delete maxwell_;
  delete vlasov_i_;
  delete vlasov_e_;
}


int VeViM::get_scheme_id() { return parameter_.kVeViM; }


real VeViM::get_dt(int species) {
  if (species == Parameter::kIon) {
    return dt_*electron_substeps_;
  }
  return dt_;
}


int VeViM::get_data_size() {
  // f_e (ncells_xv_e), f_i (ncells_xv_i)
  // E (3*ncells), B (3*ncells), B_ipol (3*ncells), maxwell_dt (1), maxwell_subcycles (1)
  return parameter_.species[parameter_.kElectron].ncells_xv +
         parameter_.species[parameter_.kIon].ncells_xv +
         9*parameter_.ncells_x + 2;
}


void VeViM::init() {

  if (parameter_.restart) {

    real scalars[3];
    restart::load_scalars(scalars, parameter_.output_directory, block_id_.my_coords);
    maxwell_->dt = scalars[0];
    maxwell_->substeps = round(scalars[1]);
    electron_substep_counter_ = round(scalars[2]);

    restart::load_field(vlasov_e_->f, 0, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(vlasov_i_->f, 1, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(maxwell_->E, 2, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(maxwell_->B, 3, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(maxwell_->B_ipol, 4, parameter_.output_directory, block_id_.my_coords);
    #pragma acc update device(vlasov_e_->f[0:parameter_.species[parameter_.kElectron].ncells_xv], \
                              vlasov_i_->f[0:parameter_.species[parameter_.kIon].ncells_xv], \
                              maxwell_->E[0:3*parameter_.ncells_x], \
                              maxwell_->B[0:3*parameter_.ncells_x], \
                              maxwell_->B_ipol[0:3*parameter_.ncells_x])
  } else {
    // load setup
    model::Model *plasma_models[2] =  {vlasov_e_, vlasov_i_};
    setup::init(plasma_models, maxwell_, parameter_);

    // calc dt
    dt_ = std::min(vlasov_e_->get_max_dt(), vlasov_i_->get_max_dt());
    if (!parameter_.relax_fdtd_cfl) dt_ = std::min(dt_, maxwell_->get_max_dt());
    dt_ = boundary_.mpi_min(dt_);
    maxwell_->dt = std::min(maxwell_->get_max_dt(), dt_ / parameter_.initial_maxwell_steps_per_fluid_step);
    maxwell_->substeps = floor(dt_ / maxwell_->dt);
    maxwell_->substeps = maxwell_->substeps - (maxwell_->substeps + 1) % 2; // force uneven

    // calc j
    convert::f_to_j(j_, vlasov_e_->f, vlasov_i_->f, parameter_);
    std::string orientation_j[3] = {"ccc", "ccc", "ccc"};
    boundary_.exchange_field(j_, parameter_.bd_cond_v, orientation_j, 3);
    ipol::yee_centered_to_face(j_, j_fc_, parameter_);
    if (parameter_.cweno_calc_j_from_flux) {
      std::cerr<<std::endl<<"Warning: Parameter::cweno_calc_j_from_flux==true is not "
                <<"available with a pure Vlasov model. Set to false for coupling."<<std::endl<<std::endl;
    }

    // maxwell init step
    std::string orientation_E[3] = {"fcc", "cfc", "ccf"};
    std::string orientation_B[3] = {"cff", "fcf", "ffc"};
    std::string orientation_EB_cc[3] = {"ccc", "ccc", "ccc"};

    maxwell_->step_E(0.5 * maxwell_->dt, j_fc_);
    boundary_.exchange_field(maxwell_->E, parameter_.bd_cond_E, orientation_E, 3);
    for (int i = 0; i < maxwell_->substeps / 2; ++i) {
      maxwell_->step_B(maxwell_->dt);
      boundary_.exchange_field(maxwell_->B, parameter_.bd_cond_B, orientation_B, 3);
      maxwell_->step_E(maxwell_->dt, j_fc_);
      boundary_.exchange_field(maxwell_->E, parameter_.bd_cond_E, orientation_E, 3);
    }
    maxwell_->step_B_ipol(maxwell_->dt);
    boundary_.exchange_field(maxwell_->B, parameter_.bd_cond_B, orientation_B, 3);
    boundary_.exchange_field(maxwell_->B_ipol, parameter_.bd_cond_B, orientation_B, 3);
  }

  dt_ = maxwell_->substeps * maxwell_->dt;
  if (parameter_.electron_subcycling) {
    electron_substeps_ = std::max((int) (vlasov_i_->get_max_dt()/dt_), 1);
    electron_substeps_ = boundary_.mpi_min(electron_substeps_);
    electron_substeps_ = electron_substeps_ - (electron_substeps_ + 1) % 2; // force uneven
  }
}


void VeViM::step() {

  ion_step_ = (electron_substep_counter_ == electron_substeps_/2);
  real dt_i = dt_*electron_substeps_;

  ipol::yee_face_to_centered(maxwell_->E, E_cc_, parameter_);
  ipol::yee_edge_to_centered(maxwell_->B_ipol, B_cc_, parameter_);

  // needed for Vlasov's first half-step v (TODO: only exchange for V-V and V-F10 neighbours)
  std::string orientation_EB_cc[3] = {"ccc", "ccc", "ccc"};
  boundary_.exchange_field(E_cc_, parameter_.bd_cond_E, orientation_EB_cc, 3);
  boundary_.exchange_field(B_cc_, parameter_.bd_cond_B, orientation_EB_cc, 3);

  // vlasov
  // half step v
  if (parameter_.vlasov_solver == "pfc") { // backsubstitution method needs no additional splitting
    vlasov_e_->step_vx(dt_ * 0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    if (ion_step_) vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    if (parameter_.dimensionality_v > 1) {
      vlasov_e_->step_vy(dt_*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (ion_step_) vlasov_i_->step_vy(dt_i*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
    if (parameter_.dimensionality_v > 2) {
      vlasov_e_->step_vz(dt_*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (ion_step_) vlasov_i_->step_vz(dt_i*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
  }
  else { // splitting for other solvers
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (ion_step_) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (ion_step_) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
    vlasov_e_->step_vx(dt_ * 0.5, *E_cc_, *B_cc_, include_boundary_cells_); // 1. x
    if (ion_step_) vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (ion_step_) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (ion_step_) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
  }


  // step x
  vlasov_e_->step_x(dt_);
  if (ion_step_) vlasov_i_->step_x(dt_i);
  if (parameter_.dimensionality_x > 1) {
    vlasov_e_->step_y(dt_);
    if (ion_step_) vlasov_i_->step_y(dt_i);
  }
  if (parameter_.dimensionality_x > 2) {
    vlasov_e_->step_z(dt_);
    if (ion_step_) vlasov_i_->step_z(dt_i);
  }


  // second half step v
  if (parameter_.vlasov_solver == "pfc") { // backsubstitution method needs no additional splitting
    vlasov_e_->step_vx(dt_ * 0.5, *E_cc_, *B_cc_);
    if (ion_step_) vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_);
    if (parameter_.dimensionality_v > 1) {
      vlasov_e_->step_vy(dt_*0.5, *E_cc_, *B_cc_);
      if (ion_step_) vlasov_i_->step_vy(dt_i*0.5, *E_cc_, *B_cc_);
    }
    if (parameter_.dimensionality_v > 2) {
      vlasov_e_->step_vz(dt_*0.5, *E_cc_, *B_cc_);
      if (ion_step_) vlasov_i_->step_vz(dt_i*0.5, *E_cc_, *B_cc_);
    }
  }
  else { // splitting for other solvers
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *E_cc_, *B_cc_);
      if (ion_step_) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_);
    }
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *E_cc_, *B_cc_);
      if (ion_step_) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_);
    }
    vlasov_e_->step_vx(dt_ * 0.5, *E_cc_, *B_cc_); // 1. x
    if (ion_step_) vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_);
    if (parameter_.dimensionality_v > 1) { // 0.5 y
      vlasov_e_->step_vy(dt_*0.5*0.5, *E_cc_, *B_cc_);
      if (ion_step_) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_);
    }
    if (parameter_.dimensionality_v > 2) { // 0.5 z
      vlasov_e_->step_vz(dt_*0.5*0.5, *E_cc_, *B_cc_);
      if (ion_step_) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_);
    }
  }

  exchange_plasma_solvers(vlasov_e_->f, vlasov_i_->f);

  // calc j
  convert::f_to_j(j_, vlasov_e_->f, vlasov_i_->f, parameter_);
  std::string orientation_j[3] = {"ccc", "ccc", "ccc"};
  boundary_.exchange_field(j_, parameter_.bd_cond_v, orientation_j, 3);
  ipol::yee_centered_to_face(j_, j_fc_, parameter_);


  electron_substep_counter_++;

  // update dt
  int substeps_maxwell_old = maxwell_->substeps; // for change of maxwell subcycles

  if (!parameter_.freeze_dt && electron_substep_counter_ == electron_substeps_) {
    electron_substep_counter_ = 0;
    dt_ = std::min(vlasov_e_->get_max_dt(), vlasov_i_->get_max_dt());
    if (!parameter_.relax_fdtd_cfl) dt_ = std::min(dt_, maxwell_->get_max_dt());
    dt_ = boundary_.mpi_min(dt_);
    maxwell_->substeps = floor(dt_ / maxwell_->dt);
    maxwell_->substeps = maxwell_->substeps - (maxwell_->substeps + 1) % 2; // force uneven
    dt_ = maxwell_->substeps * maxwell_->dt;
    if (maxwell_->substeps <= 0) {
      std::cerr << std::endl << "Error: dt_fluid < dt_maxwell not allowed, adjust "
                << "parameter::initial_maxwell_steps_per_fluid_step or check "
                << "numerical stability!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (parameter_.electron_subcycling) {
      electron_substeps_ = std::max((int) (vlasov_i_->get_max_dt()/dt_), 1);
      electron_substeps_ = boundary_.mpi_min(electron_substeps_);
      electron_substeps_ = electron_substeps_ - (electron_substeps_ + 1) % 2; // force uneven
    }
  }

  // maxwell
  std::string orientation_E[3] = {"fcc", "cfc", "ccf"};
  std::string orientation_B[3] = {"cff", "fcf", "ffc"};

  maxwell_->step_E(maxwell_->dt, j_fc_);
  boundary_.exchange_field(maxwell_->E, parameter_.bd_cond_E, orientation_E, 3);
  for (int i = 0; i < (substeps_maxwell_old + maxwell_->substeps) / 2 - 1; ++i) {
    maxwell_->step_B(maxwell_->dt);
    boundary_.exchange_field(maxwell_->B, parameter_.bd_cond_B, orientation_B, 3);
    maxwell_->step_E(maxwell_->dt, j_fc_);
    boundary_.exchange_field(maxwell_->E, parameter_.bd_cond_E, orientation_E, 3);
  }
  maxwell_->step_B_ipol(maxwell_->dt);
  boundary_.exchange_field(maxwell_->B, parameter_.bd_cond_B, orientation_B, 3);
  boundary_.exchange_field(maxwell_->B_ipol, parameter_.bd_cond_B, orientation_B, 3);
}


void VeViM::exchange_plasma_solvers(real* data_e, real* data_i) {

  real* buffer_out_e[parameter_.kNDirections],
      * buffer_out_i[parameter_.kNDirections],
      * buffer_in_e [parameter_.kNDirections],
      * buffer_in_i [parameter_.kNDirections];
  int buffer_size_e[parameter_.kNDirections],
      buffer_size_i[parameter_.kNDirections];

  std::string orientation[1] = {"ccc"};

  bool exchange_this_direction_e[6] = {0, 0, 0, 0, 0, 0};
  bool exchange_this_direction_i[6] = {0, 0, 0, 0, 0, 0};

  for (int d = 0; d < parameter_.kNDirections; ++d) {

    // prepare buffers
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViM) {
      // electrons
      exchange_this_direction_e[d] = true;
      buffer_size_e[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kElectron);
      buffer_out_e[d] = new real[buffer_size_e[d]];
      buffer_in_e[d]  = new real[buffer_size_e[d]];
      #pragma acc enter data create(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                    buffer_in_e [d:1][0:buffer_size_e[d]])
      boundary_.extract_buffer(data_e, buffer_out_e[d], d, 1, true, parameter_.kElectron);
      // ions
      if (ion_step_) {
        exchange_this_direction_i[d] = true;
        buffer_size_i[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kIon);
        buffer_out_i[d] = new real[buffer_size_i[d]];
        buffer_in_i[d]  = new real[buffer_size_i[d]];
        #pragma acc enter data create(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                      buffer_in_i [d:1][0:buffer_size_i[d]])
        boundary_.extract_buffer(data_i, buffer_out_i[d], d, 1, true, parameter_.kIon);
      }
    }
    else if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViM) {
      // electrons
      exchange_this_direction_e[d] = true;
      buffer_size_e[d] = boundary_.calc_buffer_size(d, 10);
      buffer_out_e[d] = new real[buffer_size_e[d]];
      buffer_in_e[d]  = new real[buffer_size_e[d]];
      #pragma acc enter data create(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                    buffer_in_e [d:1][0:buffer_size_e[d]])
      int buffer_size_tmp = boundary_.calc_buffer_size(d, 1, true, parameter_.kElectron);
      real* buffer_tmp = new real[buffer_size_tmp];
      #pragma acc enter data create(buffer_tmp[buffer_size_tmp])
      boundary_.extract_buffer(data_e, buffer_tmp, d, 1, true, parameter_.kElectron);
      convert::f_to_ten_moments(buffer_out_e[d], buffer_tmp,
          boundary_.calc_buffer_array_dimensions(d),
          parameter_.species[parameter_.kElectron], parameter_);
      #pragma acc exit data delete(buffer_tmp[buffer_size_tmp])
      delete[] buffer_tmp;
      // ions
      if (ion_step_) {
        exchange_this_direction_i[d] = true;
        buffer_size_i[d] = boundary_.calc_buffer_size(d, 1, true, parameter_.kIon);
        buffer_out_i[d] = new real[buffer_size_i[d]];
        buffer_in_i[d]  = new real[buffer_size_i[d]];
        #pragma acc enter data create(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                      buffer_in_i [d:1][0:buffer_size_i[d]])
        boundary_.extract_buffer(data_i, buffer_out_i[d], d, 1, true, parameter_.kIon);
      }
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

    // apply received buffers (convert if necessary)
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kVeViM) {
      // electrons
      boundary_.apply_buffer(data_e, buffer_in_e[d], d, 1, true, parameter_.kElectron);
      // ions
      if (ion_step_) {
        boundary_.apply_buffer(data_i, buffer_in_i[d], d, 1, true, parameter_.kIon);
      }
    }
    else if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViM) {
      // electrons
      int buffer_size_tmp = boundary_.calc_buffer_size(d, 1, true, parameter_.kElectron);
      real* buffer_tmp = new real[buffer_size_tmp];
      #pragma acc enter data create(buffer_tmp[0:buffer_size_tmp])
      convert::ten_moments_to_f(buffer_tmp, buffer_in_e[d],
              boundary_.calc_buffer_array_dimensions(d),
              parameter_.species[parameter_.kElectron], parameter_);
      boundary_.apply_buffer(data_e, buffer_tmp, d, 1, true, parameter_.kElectron);
      #pragma acc exit data delete(buffer_tmp[0:buffer_size_tmp])
      delete[] buffer_tmp;
      if (parameter_.ten_moment_vlasov_coupling_fit) {
        convert::delta_f_fit(data_e, buffer_in_e[d], boundary_.calc_buffer_array_dimensions(d),
            parameter_.species[parameter_.kElectron], parameter_, d);
      }
      // ions
      if (ion_step_) {
        boundary_.apply_buffer(data_i, buffer_in_i[d], d, 1, true, parameter_.kIon);
      }
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
    boundary_.apply_boundary_conditions(data_e, parameter_.bd_cond_v,
            orientation, d, 1, true, parameter_.kElectron);
    if (ion_step_) {
      boundary_.apply_boundary_conditions(data_i, parameter_.bd_cond_v,
              orientation, d, 1, true, parameter_.kIon);
    }
  }
}


void VeViM::send_model_data(real* out_buffer, int receiver_id) {
  
  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  load_balancing::send_model_data(out_buffer, get_data_size(), plasma_models, maxwell_, 2,
                                  receiver_id, boundary_, parameter_);
}


void VeViM::receive_model_data(real* in_buffer, int sender_id) {

  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  load_balancing::receive_model_data(in_buffer, get_data_size(), plasma_models, maxwell_, 2,
                                     sender_id, boundary_, parameter_);
}


void VeViM::replace_dealloc_old_models(model::Model* plasma_model[],
                  model::Model* electromagnetic_model) {
  // replace models when the scheme is changed
  delete maxwell_;
  delete vlasov_i_;
  delete vlasov_e_;

  vlasov_e_ = (model::Vlasov*) plasma_model[parameter_.kElectron];
  vlasov_i_ = (model::Vlasov*) plasma_model[parameter_.kIon];
  maxwell_ = (model::Maxwell*) electromagnetic_model;
}


void VeViM::calc_alloc_converted_model(model::Model*& model_out,
                  int species, int target_scheme_id) {
  // convert models when the scheme is changed

  switch (species) {
    case Parameter::kElectron: {
      if (target_scheme_id == parameter_.kF10eViM || target_scheme_id == parameter_.kF10eF10iM) { // to 10 mom
        model_out = new model::Fluid10(parameter_, parameter_.species[parameter_.kElectron]);
        model::Fluid10* model_cast = (model::Fluid10*) model_out;
        convert::f_to_ten_moments(model_cast->ten_moments, vlasov_e_->f,
                  parameter_.species[parameter_.kElectron], parameter_);
      } else if (target_scheme_id == parameter_.kF5eF10iM || target_scheme_id == parameter_.kF5eF5iM) { // to 5 mom
        model_out = new model::Fluid5(parameter_, parameter_.species[parameter_.kElectron]);
        model::Fluid5* model_cast = (model::Fluid5*) model_out;
        convert::f_to_five_moments(model_cast->five_moments, vlasov_e_->f,
                  parameter_.species[parameter_.kElectron], parameter_);
      }
    } break;
    case Parameter::kIon: {
      if (target_scheme_id == parameter_.kF10eViM) { // to f
        model_out = new model::Vlasov(parameter_, parameter_.species[parameter_.kIon]);
        model::Vlasov* model_cast = (model::Vlasov*) model_out;
        util::copy_array(model_cast->f, vlasov_i_->f, parameter_.species[parameter_.kIon].ncells_xv);
      } else if (target_scheme_id == parameter_.kF10eF10iM || target_scheme_id == parameter_.kF5eF10iM) { // to 10 mom
        model_out = new model::Fluid10(parameter_, parameter_.species[parameter_.kIon]);
        model::Fluid10* model_cast = (model::Fluid10*) model_out;
        convert::f_to_ten_moments(model_cast->ten_moments, vlasov_i_->f,
                  parameter_.species[parameter_.kIon], parameter_);
      } else if (target_scheme_id == parameter_.kF5eF5iM) { // to 5 mom
        model_out = new model::Fluid5(parameter_, parameter_.species[parameter_.kIon]);
        model::Fluid5* model_cast = (model::Fluid5*) model_out;
        convert::f_to_five_moments(model_cast->five_moments, vlasov_i_->f,
                  parameter_.species[parameter_.kIon], parameter_);
      }
    } break;
    case -1: { // electromagnetic model
      model_out = new model::Maxwell(parameter_);
      model::Maxwell* model_cast = (model::Maxwell* ) model_out;
      util::copy_array(model_cast->E, maxwell_->E, 3*parameter_.ncells_x);
      util::copy_array(model_cast->B, maxwell_->B, 3*parameter_.ncells_x);
      util::copy_array(model_cast->B_ipol, maxwell_->B_ipol, 3*parameter_.ncells_x);
      model_cast->dt = maxwell_->dt;
      model_cast->substeps = maxwell_->substeps;
    } break;
  }
}


void VeViM::output(real t, int output_number, bool output_vtk) {
  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  output::prepare_and_write(plasma_models, maxwell_, 2, t, output_number, output_vtk,
                            block_id_, parameter_, electron_substep_counter_);
}


int VeViM::evaluate_criterion() {
  model::Model *plasma_models[2] = {vlasov_e_, vlasov_i_};
  return criterion::evaluate_criterion(plasma_models, maxwell_, block_id_, parameter_);
}

} // namespace scheme

/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "F10eViMdeltaf.h"

#include <iostream>
#include <cstdlib> // exit
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

// ten moments fluid electrons, Vlasov ions, Maxwell

F10eViMdeltaf::F10eViMdeltaf(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter)
    : block_id_(block_id), boundary_(mpi_boundary), parameter_(parameter) {

  fluid_e_  = new model::Fluid10(parameter_, parameter_.species[parameter_.kElectron]);
  fluid_i_  = new model::Fluid10(parameter_, parameter_.species[parameter_.kIon]);
  vlasov_i_ = new model::Vlasov(parameter_, parameter_.species[parameter_.kIon]);
  maxwell_  = new model::Maxwell(parameter_);

  un_e_ = new real[3 * parameter_.ncells_x]();
  un_i_ = new real[3 * parameter_.ncells_x]();
  j_fc_ = new real[3 * parameter_.ncells_x]();
  E_cc_ = new real[3 * parameter_.ncells_x]();
  B_cc_ = new real[3 * parameter_.ncells_x]();

  Q_raw_old_i_ = new real[10 * parameter_.ncells_x]();
  Q_raw_new_i_ = new real[10 * parameter_.ncells_x]();

  ten_moments_tmp_e_ = new real[10 * parameter_.ncells_x]();
  ten_moments_tmp_i_ = new real[10 * parameter_.ncells_x]();
  total_flux_e_ = new real[10 * parameter_.ncells_x_nobd]();
  total_flux_i_ = new real[10 * parameter_.ncells_x_nobd]();
  source_e_ = new real[10 * parameter_.ncells_x_nobd]();
  source_i_ = new real[10 * parameter_.ncells_x_nobd]();

  #pragma acc enter data copyin(un_e_[0:3*parameter_.ncells_x], \
                                un_i_[0:3*parameter_.ncells_x], \
                                j_fc_[0:3*parameter_.ncells_x], \
                                E_cc_[0:3*parameter_.ncells_x], \
                                B_cc_[0:3*parameter_.ncells_x], \
                                Q_raw_old_i_[0:10*parameter_.ncells_x], \
                                Q_raw_new_i_[0:10*parameter_.ncells_x], \
                                ten_moments_tmp_e_[0:10*parameter_.ncells_x], \
                                ten_moments_tmp_i_[0:10*parameter_.ncells_x], \
                                total_flux_e_[0:10*parameter_.ncells_x_nobd], \
                                total_flux_i_[0:10*parameter_.ncells_x_nobd], \
                                source_e_[0:10*parameter_.ncells_x_nobd], \
                                source_i_[0:10*parameter_.ncells_x_nobd])
}


F10eViMdeltaf::~F10eViMdeltaf() {
  #pragma acc exit data delete(un_e_[0:3*parameter_.ncells_x], \
                               un_i_[0:3*parameter_.ncells_x], \
                               j_fc_[0:3*parameter_.ncells_x], \
                               E_cc_[0:3*parameter_.ncells_x], \
                               B_cc_[0:3*parameter_.ncells_x], \
                               Q_raw_old_i_[0:10*parameter_.ncells_x], \
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
  delete[] Q_raw_old_i_;

  delete[] B_cc_;
  delete[] E_cc_;
  delete[] j_fc_;
  delete[] un_i_;
  delete[] un_e_;

  delete maxwell_;
  delete fluid_i_;
  delete fluid_e_;
  delete vlasov_i_;
}


int F10eViMdeltaf::get_scheme_id() { return parameter_.kF10eViMdeltaf; }


real F10eViMdeltaf::get_dt(int species) {
  if (species == Parameter::kIon) {
    return dt_*electron_substeps_;
  }
  return dt_;
}


int F10eViMdeltaf::get_data_size() {
  // f_i (ncells_xv_i)
  // ten_moments_e (10*ncells), ten_moments_i (10*ncells)
  // E (3*ncells), B (3*ncells), B_ipol (3*ncells), maxwell_dt (1), maxwell_subcycles (1)
  return parameter_.species[parameter_.kIon].ncells_xv + 29*parameter_.ncells_x + 2;
}


void F10eViMdeltaf::init() {

  if (parameter_.restart) {

    real scalars[3];
    restart::load_scalars(scalars, parameter_.output_directory, block_id_.my_coords);
    maxwell_->dt = scalars[0];
    maxwell_->substeps = round(scalars[1]);
    electron_substep_counter_ = round(scalars[2]);

    restart::load_field(fluid_e_->ten_moments, 0, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(vlasov_i_->f, 1, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(fluid_i_->ten_moments, 2, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(maxwell_->E, 3, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(maxwell_->B, 4, parameter_.output_directory, block_id_.my_coords);
    restart::load_field(maxwell_->B_ipol, 5, parameter_.output_directory, block_id_.my_coords);

    #pragma acc update device(fluid_e_->ten_moments[0:10*parameter_.ncells_x], \
                              vlasov_i_->f[0:parameter_.species[parameter_.kIon].ncells_xv], \
                              fluid_i_->ten_moments[0:10*parameter_.ncells_x], \
                              maxwell_->E[0:3*parameter_.ncells_x], \
                              maxwell_->B[0:3*parameter_.ncells_x], \
                              maxwell_->B_ipol[0:3*parameter_.ncells_x])
  } else {
    // load setup
    model::Model *plasma_models[2] = {fluid_e_, vlasov_i_};
    setup::init(plasma_models, maxwell_, parameter_);
    convert::f_to_ten_moments(fluid_i_->ten_moments, vlasov_i_->f, parameter_.species[parameter_.kIon], parameter_);

    // calc dt
    real dt_vlasov = vlasov_i_->get_max_dt();
    real dt_fluid = std::min(fluid_e_->get_max_dt(), fluid_i_->get_max_dt());
    dt_ = std::min(dt_vlasov, dt_fluid);
    if (!parameter_.relax_fdtd_cfl) dt_ = std::min(dt_, maxwell_->get_max_dt());
    dt_ = boundary_.mpi_min(dt_);
    maxwell_->dt = std::min(maxwell_->get_max_dt(), dt_ / parameter_.initial_maxwell_steps_per_fluid_step);
    maxwell_->substeps = floor(dt_ / maxwell_->dt);
    maxwell_->substeps = maxwell_->substeps - (maxwell_->substeps + 1) % 2; // force uneven

    // calc j
    real *j = new real[3 * parameter_.ncells_x]();
    #pragma acc enter data copyin(j[0:3*parameter_.ncells_x])
    convert::un_to_j(j, &fluid_e_->ten_moments[parameter_.ncells_x],
                     &fluid_i_->ten_moments[parameter_.ncells_x],
                     parameter_); // ten_moments(1:3) = n u
    std::string orientation_j[3] = {"ccc", "ccc", "ccc"};
    boundary_.exchange_field(j, parameter_.bd_cond_v, orientation_j, 3);
    ipol::yee_centered_to_face(j, j_fc_, parameter_);
    #pragma acc exit data delete(j[0:3*parameter_.ncells_x])
    delete[] j;

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

  util::copy_array(ten_moments_tmp_e_, fluid_e_->ten_moments, 10*parameter_.ncells_x);
  util::copy_array(ten_moments_tmp_i_, fluid_i_->ten_moments, 10*parameter_.ncells_x);

  fluid_i_->calc_un_from_flux(un_i_, true);

  // calc heat flux
  std::string orientation_Q_h[10] = {"ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc"};
  bool exchange_in_this_direction_i[6] = {0, 0, 0, 0, 0, 0};
  for (int d = 0; d < parameter_.kNDirections; ++d) {
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf ||
        block_id_.neighbour_scheme_ids[d] == parameter_.kVeViMdeltaf) {
      exchange_in_this_direction_i[d] = true;
    }
  }
  convert::f_to_raw_heatflux(Q_raw_new_i_, vlasov_i_->f, parameter_.species[parameter_.kIon],
                             exchange_in_this_direction_i, parameter_);
  boundary_.exchange_field(Q_raw_new_i_, parameter_.bd_cond_Q_h, orientation_Q_h, exchange_in_this_direction_i, 10);


  dt_ = maxwell_->substeps * maxwell_->dt;
  if (parameter_.electron_subcycling) {
    electron_substeps_ = std::max((int) (std::min(fluid_i_->get_max_dt(), vlasov_i_->get_max_dt())/dt_), 1);
    electron_substeps_ = boundary_.mpi_min(electron_substeps_);
    electron_substeps_ = electron_substeps_ - (electron_substeps_ + 1) % 2; // force uneven
  }
}


void F10eViMdeltaf::step() {

  ion_step_ = (electron_substep_counter_ == electron_substeps_/2);
  real dt_i = dt_*electron_substeps_;

  ipol::yee_face_to_centered(maxwell_->E, E_cc_, parameter_);
  ipol::yee_edge_to_centered(maxwell_->B_ipol, B_cc_, parameter_);

  // needed for Vlasov's first half-step v (TODO only exchange for V-V and V-F10 neighbours)
  std::string orientation_EB_cc[3] = {"ccc", "ccc", "ccc"};
  boundary_.exchange_field(E_cc_, parameter_.bd_cond_E, orientation_EB_cc, 3);
  boundary_.exchange_field(B_cc_, parameter_.bd_cond_B, orientation_EB_cc, 3);


  if (ion_step_) {
    util::copy_array(Q_raw_old_i_, Q_raw_new_i_, 10*parameter_.ncells_x);

    // vlasov ions half step v
    if (parameter_.vlasov_solver == "pfc") { // backsubstitution method needs no additional splitting
      vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (parameter_.dimensionality_v > 1) vlasov_i_->step_vy(dt_i*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
      if (parameter_.dimensionality_v > 2) vlasov_i_->step_vz(dt_i*0.5, *E_cc_, *B_cc_, include_boundary_cells_);
    }
    else { // splitting for other solvers
      if (parameter_.dimensionality_v > 2) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_); // 0.5 z
      if (parameter_.dimensionality_v > 1) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);  // 0.5 y
      vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_, include_boundary_cells_);                                         // 1. x
      if (parameter_.dimensionality_v > 1) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_);  // 0.5 y
      if (parameter_.dimensionality_v > 2) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_, include_boundary_cells_); // 0.5 z
    }

    // vlasov ions step x
    vlasov_i_->step_x(dt_i);
    if (parameter_.dimensionality_x > 1) vlasov_i_->step_y(dt_i);
    if (parameter_.dimensionality_x > 2) vlasov_i_->step_z(dt_i);

    // vlasov ions second half step v
    if (parameter_.vlasov_solver == "pfc") { // backsubstitution method needs no additional splitting
      vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_);
      if (parameter_.dimensionality_v > 1) vlasov_i_->step_vy(dt_i*0.5, *E_cc_, *B_cc_);
      if (parameter_.dimensionality_v > 2) vlasov_i_->step_vz(dt_i*0.5, *E_cc_, *B_cc_);
    }
    else { // splitting for other solvers
      if (parameter_.dimensionality_v > 2) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_); // 0.5 z
      if (parameter_.dimensionality_v > 1) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_);  // 0.5 y
      vlasov_i_->step_vx(dt_i * 0.5, *E_cc_, *B_cc_);                                         // 1. x
      if (parameter_.dimensionality_v > 1) vlasov_i_->step_vy(dt_i*0.5*0.5, *E_cc_, *B_cc_);  // 0.5 y
      if (parameter_.dimensionality_v > 2) vlasov_i_->step_vz(dt_i*0.5*0.5, *E_cc_, *B_cc_); // 0.5 z
    }


    // calc heat flux
    std::string orientation_Q_h[10] = {"ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc", "ccc"};
    bool exchange_in_this_direction_i[6] = {0, 0, 0, 0, 0, 0};
    for (int d = 0; d < parameter_.kNDirections; ++d) {
      if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf ||
          block_id_.neighbour_scheme_ids[d] == parameter_.kVeViMdeltaf) {
        exchange_in_this_direction_i[d] = true;
      }
    }
    convert::f_to_raw_heatflux(Q_raw_new_i_, vlasov_i_->f, parameter_.species[parameter_.kIon],
                               exchange_in_this_direction_i, parameter_);
    boundary_.exchange_field(Q_raw_new_i_, parameter_.bd_cond_Q_h, orientation_Q_h, exchange_in_this_direction_i, 10);
  } // end if ion_step_

  // fluid rk step 1
  fluid_e_->calc_single_step_rk3(ten_moments_tmp_e_, total_flux_e_, source_e_, E_cc_, B_cc_);
  heat_flux_closure(fluid_e_, ten_moments_tmp_e_, source_e_);
  fluid_e_->apply_single_step_rk3(dt_, ten_moments_tmp_e_, total_flux_e_, source_e_, 1., 0., 1.);
  if (ion_step_) {
    fluid_i_->calc_single_step_rk3(ten_moments_tmp_i_, total_flux_i_, source_i_, Q_raw_old_i_, E_cc_, B_cc_);
    fluid_i_->apply_single_step_rk3(dt_i, ten_moments_tmp_i_, total_flux_i_, source_i_, 1., 0., 1.);
  }
  exchange_plasma_solvers(vlasov_i_->f, ten_moments_tmp_e_, ten_moments_tmp_i_, runge_kutta_only_);

  // fluid rk step 2
  fluid_e_->calc_single_step_rk3(ten_moments_tmp_e_, total_flux_e_, source_e_, E_cc_, B_cc_);
  heat_flux_closure(fluid_e_, ten_moments_tmp_e_, source_e_);
  fluid_e_->apply_single_step_rk3(dt_, ten_moments_tmp_e_, total_flux_e_, source_e_, .75, .25, .25);
  if (ion_step_) {
    fluid_i_->calc_single_step_rk3(ten_moments_tmp_i_, total_flux_i_, source_i_, Q_raw_new_i_, E_cc_, B_cc_);
    fluid_i_->apply_single_step_rk3(dt_i, ten_moments_tmp_i_, total_flux_i_, source_i_, .75, .25, .25);
  }
  exchange_plasma_solvers(vlasov_i_->f, ten_moments_tmp_e_, ten_moments_tmp_i_, runge_kutta_only_);

  // interpolate heat flux to intermediate time: Q_old = 0.5*Q_old + 0.5*Q_new
  if (ion_step_) {
    util::add_arrays(Q_raw_old_i_, Q_raw_old_i_, Q_raw_new_i_, 0.5, 0.5, 10*parameter_.ncells_x);
  }

  // fluid rk step 3
  fluid_e_->calc_single_step_rk3(ten_moments_tmp_e_, total_flux_e_, source_e_, E_cc_, B_cc_);
  heat_flux_closure(fluid_e_, ten_moments_tmp_e_, source_e_);
  fluid_e_->apply_single_step_rk3(dt_, ten_moments_tmp_e_, total_flux_e_, source_e_, 1./3., 2./3., 2./3.);
  if (ion_step_) {
    fluid_i_->calc_single_step_rk3(ten_moments_tmp_i_, total_flux_i_, source_i_, Q_raw_old_i_, E_cc_, B_cc_);
    fluid_i_->apply_single_step_rk3(dt_i, ten_moments_tmp_i_, total_flux_i_, source_i_, 1./3., 2./3., 2./3.);

    convert::delta_f_correction(vlasov_i_->f, ten_moments_tmp_i_,
                                parameter_.species[parameter_.kIon], parameter_);
  }

  exchange_plasma_solvers(vlasov_i_->f, ten_moments_tmp_e_, ten_moments_tmp_i_);

  util::copy_array(fluid_e_->ten_moments, ten_moments_tmp_e_, 10*parameter_.ncells_x);
  if (ion_step_) {
    util::copy_array(fluid_i_->ten_moments, ten_moments_tmp_i_, 10*parameter_.ncells_x);
  }


  // calc j
  if (parameter_.cweno_calc_j_from_flux) {
    fluid_e_->calc_un_from_flux(un_e_);
    if (ion_step_) {
      fluid_i_->calc_un_from_flux(un_i_, true);
    }
    convert::un_to_j(j_fc_, un_e_, un_i_,
                     &fluid_e_->ten_moments[parameter_.ncells_x],
                     &fluid_i_->ten_moments[parameter_.ncells_x], parameter_, block_id_);
  } else { // workaround until j from flux is implemented for Vlasov
    real *j = new real[3 * parameter_.ncells_x]();
    #pragma acc enter data copyin(j[0:3*parameter_.ncells_x])
    convert::un_to_j(j, &fluid_e_->ten_moments[parameter_.ncells_x],
                     &fluid_i_->ten_moments[parameter_.ncells_x],
                     parameter_); // ten_moments(1:3) = n u
    std::string orientation_j[3] = {"ccc", "ccc", "ccc"};
    boundary_.exchange_field(j, parameter_.bd_cond_v, orientation_j, 3);
    ipol::yee_centered_to_face(j, j_fc_, parameter_);
    #pragma acc exit data delete(j[0:3*parameter_.ncells_x])
    delete[] j;
  }

  electron_substep_counter_++;

  // update dt
  int substeps_maxwell_old = maxwell_->substeps; // for change of maxwell subcycles

  if (!parameter_.freeze_dt && electron_substep_counter_ == electron_substeps_) {
    electron_substep_counter_ = 0;
    real dt_vlasov = vlasov_i_->get_max_dt();
    real dt_fluid = std::min(fluid_e_->get_max_dt(), fluid_i_->get_max_dt());
    dt_ = std::min(dt_vlasov, dt_fluid);
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
      electron_substeps_ = std::max((int) (std::min(fluid_i_->get_max_dt(), vlasov_i_->get_max_dt())/dt_), 1);
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


void F10eViMdeltaf::heat_flux_closure(model::Fluid10 *fluid, const real *ten_moments_tmp, real *source) {
  if (parameter_.ten_moment_closure == "gradient" ||
      parameter_.ten_moment_closure == "gradient_P" ||
      parameter_.ten_moment_closure == "gradient_T" ||
      parameter_.ten_moment_closure == "gradient_sym") {
    real *P_m_tmp = new real[6 * parameter_.ncells_x];
    #pragma acc enter data create(P_m_tmp[0:6*parameter_.ncells_x])
    std::string orientation_P[6] = {"ccc", "ccc", "ccc", "ccc", "ccc", "ccc"};

    bool exchange_in_this_direction[6] = {0, 0, 0, 0, 0, 0};
    for (int d = 0; d < parameter_.kNDirections; ++d) {
      if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eF10iM ||
         block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf) {
        exchange_in_this_direction[d] = true;
      }
    }

    for (int i = 0; i < parameter_.gradient_closure_subcycles; ++i) { // subcycling of laplace P
      fluid->gradient_heat_flux_closure_cycle(i, dt_, ten_moments_tmp, source, P_m_tmp);

      if (i < parameter_.gradient_closure_subcycles - 1) {
        boundary_.exchange_field(P_m_tmp, parameter_.bd_cond_P, orientation_P,
                exchange_in_this_direction, 6);
      }
    }

    #pragma acc exit data delete(P_m_tmp[0:6*parameter_.ncells_x])
    delete[] P_m_tmp;
  }

  else if (parameter_.ten_moment_closure == "isotropization" ||
           parameter_.ten_moment_closure == "scalar_k") { // scalar_k is a legacy name
    fluid->isotropization_heat_flux_closure(ten_moments_tmp, source);
  }

  else if (parameter_.ten_moment_closure != "none") {
    std::cerr << std::endl
              << "Error: The chosen heat flux closure does not exist. "
              << "Please check the value in parameter::ten_moment_closure."
              << std::endl;
  }
}


void F10eViMdeltaf::exchange_plasma_solvers(real* data_vlasov_i,
          real* data_fluid_e, real* data_fluid_i, bool runge_kutta_only) {

  real* buffer_out_e[parameter_.kNDirections],
      * buffer_out_i[parameter_.kNDirections],
      * buffer_in_e [parameter_.kNDirections],
      * buffer_in_i [parameter_.kNDirections];
  int buffer_size_e[parameter_.kNDirections],
      buffer_size_i[parameter_.kNDirections];

  // fluid exchange
  std::string orientation[10] = {"ccc","ccc","ccc","ccc","ccc","ccc","ccc","ccc","ccc","ccc"};

  bool exchange_this_direction_e[6] = {0, 0, 0, 0, 0, 0};
  bool exchange_this_direction_i[6] = {0, 0, 0, 0, 0, 0};

  for (int d = 0; d < parameter_.kNDirections; ++d) {

    // prepare buffers
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf ||
        block_id_.neighbour_scheme_ids[d] == parameter_.kVeViMdeltaf ||
        block_id_.neighbour_scheme_ids[d] == parameter_.kF10eF10iM) {
      // electrons
      exchange_this_direction_e[d] = true;
      buffer_size_e[d] = boundary_.calc_buffer_size(d, 10);
      buffer_out_e[d] = new real[buffer_size_e[d]];
      buffer_in_e[d]  = new real[buffer_size_e[d]];
      #pragma acc enter data create(buffer_out_e[d:1][0:buffer_size_e[d]],\
                                    buffer_in_e [d:1][0:buffer_size_e[d]])
      boundary_.extract_buffer(data_fluid_e, buffer_out_e[d], d, 10);
      // ions
      if (ion_step_) {
        exchange_this_direction_i[d] = true;
        buffer_size_i[d] = boundary_.calc_buffer_size(d, 10);
        buffer_out_i[d] = new real[buffer_size_i[d]];
        buffer_in_i[d]  = new real[buffer_size_i[d]];
        #pragma acc enter data create(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                      buffer_in_i [d:1][0:buffer_size_i[d]])
        boundary_.extract_buffer(data_fluid_i, buffer_out_i[d], d, 10);
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
    if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf ||
        block_id_.neighbour_scheme_ids[d] == parameter_.kVeViMdeltaf ||
        block_id_.neighbour_scheme_ids[d] == parameter_.kF10eF10iM) {
      // electrons
      boundary_.apply_buffer(data_fluid_e, buffer_in_e[d], d, 10);

      // ions
      if (ion_step_) {
        boundary_.apply_buffer(data_fluid_i, buffer_in_i[d], d, 10);

        if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eF10iM && !runge_kutta_only) {
          int buffer_size_tmp = boundary_.calc_buffer_size(d, 1, true, parameter_.kIon);
          real* buffer_tmp = new real[buffer_size_tmp];
          #pragma acc enter data create(buffer_tmp[0:buffer_size_tmp])
          convert::ten_moments_to_f(buffer_tmp, buffer_in_i[d],
                  boundary_.calc_buffer_array_dimensions(d),
                  parameter_.species[parameter_.kIon], parameter_);
          boundary_.apply_buffer(data_vlasov_i, buffer_tmp, d, 1, true, parameter_.kIon);
          #pragma acc exit data delete(buffer_tmp[0:buffer_size_tmp])
          delete[] buffer_tmp;
          if (parameter_.ten_moment_vlasov_coupling_fit) {
            convert::delta_f_fit(data_vlasov_i, buffer_in_i[d], boundary_.calc_buffer_array_dimensions(d),
                  parameter_.species[parameter_.kIon], parameter_, d);
          }
        }
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
    boundary_.apply_boundary_conditions(data_fluid_e, parameter_.bd_cond_ten_moments,
            orientation, d, 10);
    if (ion_step_) {
      boundary_.apply_boundary_conditions(data_fluid_i, parameter_.bd_cond_ten_moments,
              orientation, d, 10);
    }
  }

  // Vlasov exchange
  if (!runge_kutta_only && ion_step_) {
    std::string orientation[1] = {"ccc"};

    bool exchange_this_direction[6] = {0, 0, 0, 0, 0, 0};

    for (int d = 0; d < parameter_.kNDirections; ++d) {

      // prepare buffers
      if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf ||
          block_id_.neighbour_scheme_ids[d] == parameter_.kVeViMdeltaf) {
        // ions
        exchange_this_direction[d] = true;
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
      if (exchange_this_direction[d])
        boundary_.exchange_buffer(buffer_out_i[d], buffer_in_i[d], buffer_size_i[d], d);
    }
    boundary_.wait(exchange_this_direction);

    for (int d = 0; d < parameter_.kNDirections; ++d) {

      // apply received buffers (convert if necessary)
      if (block_id_.neighbour_scheme_ids[d] == parameter_.kF10eViMdeltaf ||
          block_id_.neighbour_scheme_ids[d] == parameter_.kVeViMdeltaf) {
        // ions
        boundary_.apply_buffer(data_vlasov_i, buffer_in_i[d], d, 1, true, parameter_.kIon);
      }

      // clean up
      if (exchange_this_direction[d]) {
        #pragma acc exit data delete(buffer_out_i[d:1][0:buffer_size_i[d]],\
                                     buffer_in_i [d:1][0:buffer_size_i[d]])
        delete[] buffer_out_i[d];
        delete[] buffer_in_i [d];
      }

      boundary_.apply_boundary_conditions(data_vlasov_i, parameter_.bd_cond_v,
              orientation, d, 1, true, parameter_.kIon);
    }
  } // end if(!runge_kutta_only && ion_step_)
}


void F10eViMdeltaf::send_model_data(real* out_buffer, int receiver_id) {
  
  model::Model *plasma_models[3] = {fluid_e_, vlasov_i_, fluid_i_};
  load_balancing::send_model_data(out_buffer, get_data_size(), plasma_models, maxwell_, 3,
                                  receiver_id, boundary_, parameter_);
}


void F10eViMdeltaf::receive_model_data(real* in_buffer, int sender_id) {

  model::Model *plasma_models[3] = {fluid_e_, vlasov_i_, fluid_i_};
  load_balancing::receive_model_data(in_buffer, get_data_size(), plasma_models, maxwell_, 3,
                                     sender_id, boundary_, parameter_);
  util::copy_array(ten_moments_tmp_e_, fluid_e_->ten_moments, 10*parameter_.ncells_x);
  util::copy_array(ten_moments_tmp_i_, fluid_i_->ten_moments, 10*parameter_.ncells_x);
}


void F10eViMdeltaf::replace_dealloc_old_models(model::Model* plasma_model[],
                  model::Model* electromagnetic_model) {

  // TODO: adjust to contain fluid schemes as well for delta f,
  //       also check implications for /physics/plasma/criteria/select_scheme.cc

  // replace models when the scheme is changed
  delete maxwell_;
  delete vlasov_i_;
  delete fluid_e_;

  fluid_e_ = (model::Fluid10*) plasma_model[parameter_.kElectron];
  vlasov_i_ = (model::Vlasov*) plasma_model[parameter_.kIon];
  maxwell_ = (model::Maxwell*) electromagnetic_model;

  util::copy_array(ten_moments_tmp_e_, fluid_e_->ten_moments, 10*parameter_.ncells_x);
  convert::f_to_ten_moments(fluid_i_->ten_moments, vlasov_i_->f, parameter_.species[parameter_.kIon], parameter_);
  util::copy_array(ten_moments_tmp_i_, fluid_i_->ten_moments, 10*parameter_.ncells_x);
}


void F10eViMdeltaf::calc_alloc_converted_model(model::Model*& model_out,
                  int species, int target_scheme_id) {

  // TODO: adjust to contain fluid schemes as well for delta f,
  //       also check implications for /physics/plasma/criteria/select_scheme.cc

  // convert models when the scheme is changed

  switch (species) {
    case Parameter::kElectron: {
      if (target_scheme_id == parameter_.kVeViMdeltaf) { // to f
        model_out = new model::Vlasov(parameter_, parameter_.species[parameter_.kElectron]);
        model::Vlasov* model_cast = (model::Vlasov*) model_out;
        convert::ten_moments_to_f(model_cast->f, fluid_e_->ten_moments,
                  parameter_.species[parameter_.kElectron], parameter_);
      } else if (target_scheme_id == parameter_.kF10eF10iM) { // to 10 mom
        model_out = new model::Fluid10(parameter_, parameter_.species[Parameter::kElectron]);
        model::Fluid10* model_cast = (model::Fluid10*) model_out;
        util::copy_array(model_cast->ten_moments, fluid_e_->ten_moments, 10*parameter_.ncells_x);
      } else if (target_scheme_id == parameter_.kF5eF10iM || target_scheme_id == parameter_.kF5eF5iM) { // to 5 mom
        model_out = new model::Fluid5(parameter_, parameter_.species[Parameter::kElectron]);
        model::Fluid5* model_cast = (model::Fluid5*) model_out;
        convert::ten_moments_to_five_moments(model_cast->five_moments, fluid_e_->ten_moments, parameter_);
      }
    } break;
    case Parameter::kIon: {
      if (target_scheme_id == parameter_.kVeViMdeltaf) { // to f
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


void F10eViMdeltaf::output(real t, int output_number, bool output_vtk) {
  model::Model *plasma_models[3] = {fluid_e_, vlasov_i_, fluid_i_};
  output::prepare_and_write(plasma_models, maxwell_, 3, t, output_number, output_vtk,
                            block_id_, parameter_, electron_substep_counter_);
}


int F10eViMdeltaf::evaluate_criterion() {
  model::Model *plasma_models[2] = {fluid_e_, fluid_i_};
  return criterion::evaluate_criterion(plasma_models, maxwell_, block_id_, parameter_);
}

} // namespace scheme

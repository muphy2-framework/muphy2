/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "load_balancing.h"

#include <iostream>
#include <math.h>

#include "framework/parameter.h"
#include "framework/utilities.h"
#include "framework/mpi_boundary.h"
#include "framework/block.h"
#include "framework/model.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/models/fluid5.h"
#include "physics/plasma/models/mhd.h"
#include "physics/plasma/models/maxwell.h"
#include "physics/plasma/models/poisson.h"
#include "physics/plasma/models/ohm.h"


namespace load_balancing {

void send_model_data(real* out_buffer, int buffer_size, model::Model* plasma_model[],
                model::Model* electromagnetic_model, int nplasma_models, int receiver_id,
                Mpi_boundary& boundary, const Parameter& parameter) {

  int array_position = 0;

  // pack data into single array for sending,
  // copy from gpu to cpu (if applicable)

  // this function makes the following assumptions:
  // (i) only one Vlasov model per species,
  // (ii) in plasma_model[]: first most expensive model for each species,
  //      then second most expensive model for each species, and so on

  for (int i = 0; i < nplasma_models; ++i) {

    switch (plasma_model[i]->get_model_id()) {
    case Parameter::kVlasov: {
      if (i >= parameter.nspecies) {
        std::cerr<<std::endl<<"Error in load_balancing.cc: Please check order of passed models."<<std::endl;
        std::exit(EXIT_FAILURE);
      }
      model::Vlasov* model_cast = (model::Vlasov*) plasma_model[i];
      #pragma acc update host(model_cast->f[0:parameter.species[i].ncells_xv])
      util::copy_array_no_gpu(&out_buffer[array_position], model_cast->f, parameter.species[i].ncells_xv);
      array_position += parameter.species[i].ncells_xv;
    } break;
    case Parameter::kFluid5: {
      model::Fluid5* model_cast = (model::Fluid5*) plasma_model[i];
      #pragma acc update host(model_cast->five_moments[0:5*parameter.ncells_x])
      util::copy_array_no_gpu(&out_buffer[array_position], model_cast->five_moments, 5*parameter.ncells_x);
      array_position += 5*parameter.ncells_x;
    } break;
    case Parameter::kFluid10: {
      model::Fluid10* model_cast = (model::Fluid10*) plasma_model[i];
      #pragma acc update host(model_cast->ten_moments[0:10*parameter.ncells_x])
      util::copy_array_no_gpu(&out_buffer[array_position], model_cast->ten_moments, 10*parameter.ncells_x);
      array_position += 10*parameter.ncells_x;
    } break;
    case Parameter::kMHD1Temperature: {
      model::MHD1Temperature* model_cast = (model::MHD1Temperature*) plasma_model[i];
      util::copy_array_no_gpu(&out_buffer[array_position], model_cast->five_moments, 5*parameter.ncells_x);
      array_position += 5*parameter.ncells_x;
    } break;
      // other models
    }
  }

  switch (electromagnetic_model->get_model_id()) {
  case Parameter::kMaxwell: {
    model::Maxwell* model_cast = (model::Maxwell*) electromagnetic_model;
    #pragma acc update host(model_cast->E[0:3*parameter.ncells_x],\
                            model_cast->B[0:3*parameter.ncells_x],\
                            model_cast->B_ipol[0:3*parameter.ncells_x])
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->E, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->B, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->B_ipol, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    out_buffer[array_position] = model_cast->substeps;
    array_position += 1;
    out_buffer[array_position] = model_cast->dt;
    array_position += 1;
  } break;
  case Parameter::kPoisson: {
    model::Poisson* model_cast = (model::Poisson*) electromagnetic_model;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->E, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->B, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->phi, parameter.ncells_x);
    array_position += parameter.ncells_x;
  } break;
  case Parameter::kOhm: {
    model::Ohm* model_cast = (model::Ohm*) electromagnetic_model;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->E, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->B, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(&out_buffer[array_position], model_cast->J, 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
  } break;
    // other models
  }

  boundary.send_data(out_buffer, buffer_size, receiver_id);

} // end send_model_data()


void receive_model_data(real* in_buffer, int buffer_size, model::Model* plasma_model[],
                model::Model* electromagnetic_model, int nplasma_models, int sender_id,
                Mpi_boundary& boundary, const Parameter& parameter) {

  int array_position = 0;

  boundary.receive_data(in_buffer, buffer_size, sender_id);

  // unpack data into the models' arrays
  // copy from cpu to gpu (if applicable)

  // this function makes the following assumptions:
  // (i) only one Vlasov model per species,
  // (ii) in plasma_model[]: first most expensive model for each species,
  //      then second most expensive model for each species, and so on

  for (int i = 0; i < nplasma_models; ++i) {

    switch (plasma_model[i]->get_model_id()) {
    case Parameter::kVlasov: {
      if (i >= parameter.nspecies) {
        std::cerr<<std::endl<<"Error in load_balancing.cc: Please check order of passed models."<<std::endl;
        std::exit(EXIT_FAILURE);
      }
      model::Vlasov* model_cast = (model::Vlasov*) plasma_model[i];
      util::copy_array_no_gpu(model_cast->f, &in_buffer[array_position], parameter.species[i].ncells_xv);
      #pragma acc update device(model_cast->f[0:parameter.species[i].ncells_xv])
      array_position += parameter.species[i].ncells_xv;
    } break;
    case Parameter::kFluid5: {
      model::Fluid5* model_cast = (model::Fluid5*) plasma_model[i];
      util::copy_array_no_gpu(model_cast->five_moments, &in_buffer[array_position], 5*parameter.ncells_x);
      #pragma acc update device(model_cast->five_moments[0:5*parameter.ncells_x])
      array_position += 5*parameter.ncells_x;
    } break;
    case Parameter::kFluid10: {
      model::Fluid10* model_cast = (model::Fluid10*) plasma_model[i];
      util::copy_array_no_gpu(model_cast->ten_moments, &in_buffer[array_position], 10*parameter.ncells_x);
      #pragma acc update device(model_cast->ten_moments[0:10*parameter.ncells_x])
      array_position += 10*parameter.ncells_x;
    } break;
    case Parameter::kMHD1Temperature: {
      model::MHD1Temperature* model_cast = (model::MHD1Temperature*) plasma_model[i];
      util::copy_array_no_gpu(model_cast->five_moments, &in_buffer[array_position], 5*parameter.ncells_x);
      array_position += 5*parameter.ncells_x;
    } break;
      // other models
    }
  }

  switch (electromagnetic_model->get_model_id()) {
  case Parameter::kMaxwell: {
    model::Maxwell* model_cast = (model::Maxwell*) electromagnetic_model;
    util::copy_array_no_gpu(model_cast->E, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(model_cast->B, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(model_cast->B_ipol, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    model_cast->substeps = round(in_buffer[array_position]);
    array_position += 1;
    model_cast->dt = in_buffer[array_position];
    array_position += 1;
    #pragma acc update device(model_cast->E[0:3*parameter.ncells_x],\
                              model_cast->B[0:3*parameter.ncells_x],\
                              model_cast->B_ipol[0:3*parameter.ncells_x])
  } break;
  case Parameter::kPoisson: {
    model::Poisson* model_cast = (model::Poisson*) electromagnetic_model;
    util::copy_array_no_gpu(model_cast->E, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(model_cast->B, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(model_cast->phi, &in_buffer[array_position], parameter.ncells_x);
    array_position += parameter.ncells_x;
  } break;
  case Parameter::kOhm: {
    model::Ohm* model_cast = (model::Ohm*) electromagnetic_model;
    util::copy_array_no_gpu(model_cast->E, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(model_cast->B, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
    util::copy_array_no_gpu(model_cast->J, &in_buffer[array_position], 3*parameter.ncells_x);
    array_position += 3*parameter.ncells_x;
  } break;
    // other models
  }

} // end receive_model_data()


} // namespace load_balancing


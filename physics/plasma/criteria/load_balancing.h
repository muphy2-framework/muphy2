/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef LOAD_BALANCING_H_
#define LOAD_BALANCING_H_

#include "framework/definitions.h"

struct Block_id;
struct Parameter;
class Mpi_boundary;

namespace model{
class Model;
}

namespace load_balancing {

void send_model_data(real* out_buffer, int buffer_size, model::Model* plasma_model[],
                model::Model* electromagnetic_model, int nplasma_models, int receiver_id,
                Mpi_boundary& boundary, const Parameter& parameter);

void receive_model_data(real* in_buffer, int buffer_size, model::Model* plasma_model[],
                model::Model* electromagnetic_model, int nplasma_models, int sender_id,
                Mpi_boundary& boundary, const Parameter& parameter);

}

#endif // LOAD_BALANCING_H_

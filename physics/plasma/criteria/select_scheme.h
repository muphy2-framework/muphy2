/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef SELECT_SCHEME_H_
#define SELECT_SCHEME_H_

#include "framework/definitions.h"

struct Block_id;
struct Parameter;
class Mpi_boundary;

namespace scheme{
class Scheme;
}

namespace model{
class Model;
}

namespace criterion {

int evaluate_criterion(model::Model* plasma_model[], model::Model* electromagnetic_model,
        const Block_id& block_id, const Parameter& parameter);

scheme::Scheme* allocate_scheme(int new_scheme_id, const Block_id& block_id,
        Mpi_boundary& boundary, const Parameter& parameter);

void convert_scheme(scheme::Scheme* current_scheme, scheme::Scheme* new_scheme,
        const Parameter& parameter);
}

#endif // SELECT_SCHEME_H_

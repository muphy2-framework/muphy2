/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef CRITERION_J_H_
#define CRITERION_J_H_

#include "framework/criterion.h"

// criterion based on current density j

struct Parameter;

namespace model{
class Model;
}

namespace criterion {

class Criterion_j : public Criterion {
public:
  Criterion_j(int nused_schemes, int* scheme_hierarchy, real* thresholds);
  ~Criterion_j();
  int get_criterion_id();
  int evaluate(model::Model* plasma_model[], model::Model* electromagnetic_model, const Parameter& p);
private:
  int nused_schemes_;
  int* scheme_hierarchy_;
  real* thresholds_;
};

} // namespace criterion

#endif // CRITERION_J_H_

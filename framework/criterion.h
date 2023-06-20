/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef CRITERION_H_
#define CRITERION_H_

#include "framework/definitions.h"

// criterion base class

namespace criterion {

class Criterion {
public:
  virtual ~Criterion() { };
  virtual int get_criterion_id() = 0;
};

} // namespace criterion

#endif // CRITERION_H_

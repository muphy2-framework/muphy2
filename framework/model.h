/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "framework/definitions.h"

namespace model {

  // model base class for all possible models, e.g. 5 moment, Vlasov, Maxwell, ...

class Model {
public:
  virtual ~Model() { };
  virtual int get_model_id() = 0;
  virtual real get_max_dt() = 0;
};

} // namespace model

#endif // MODEL_H_

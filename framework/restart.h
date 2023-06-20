/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef RESTART_H_
#define RESTART_H_

#include <string>
#include "framework/definitions.h"

namespace restart {

  int load_scheme_id(const std::string output_directory, const int coords[]);
  real load_physical_time(const std::string output_directory, const int coords[]);
  void load_scalars(real* data, const std::string output_directory, const int coords[]);
  void load_field(real* data, int field_number, const std::string output_directory, const int coords[]);

}

#endif // RESTART_H_

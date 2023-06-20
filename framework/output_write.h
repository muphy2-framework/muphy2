/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef OUTPUT_WRITE_H_
#define OUTPUT_WRITE_H_

#include <string>
#include "framework/definitions.h"

struct Parameter;

namespace output {

  void write_vti_file(const real* fields[], const int nfields, const std::string field_names[],
            const int output_number, const std::string output_directory,
            const real xb_loc[], const real dx[], const int coords[], const int res_x[], const int bd[]);

  void write_pvti_file(const int nfields, const std::string field_names[], const real physical_time,
            const int output_number, const std::string output_directory,
            const real xb[], const real dx[], const int res_x[], const int res_x_total[]);

  void write_text_file(const real fields[], const int nfields, const std::string field_names[],
            const int output_number, const std::string output_directory, const int coords[]);

  void write_restart_file(const int scheme_id, const real physical_time,
            const real scalars[], const int nscalars,
            const real* fields[], const int field_sizes[], const int nfields,
            const std::string output_directory, const int coords[],
            const int output_number = -1);

  void write_parameter_file(const Parameter &p);
}

#endif // OUTPUT_WRITE_H_

/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "restart.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "framework/definitions.h"

namespace restart {

  int load_scheme_id(const std::string output_directory, const int coords[]) {

    std::string filepath = output_directory + "/restart/restart" +
      "_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) +
      "_" + std::to_string(coords[2]) + ".bin";

    int scheme_id = -1;

    std::ifstream is(filepath.c_str(), std::ios::binary | std::ios::in);
    is.seekg(0);
    is.read(reinterpret_cast<char*>(&scheme_id), sizeof(int));
    is.close();

    return scheme_id;
  }


  real load_physical_time(const std::string output_directory, const int coords[]) {

    std::string filepath = output_directory + "/restart/restart" +
      "_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) +
      "_" + std::to_string(coords[2]) + ".bin";

    real physical_time = 0.;

    std::ifstream is(filepath.c_str(), std::ios::binary | std::ios::in);
    is.seekg(sizeof(int));
    is.read(reinterpret_cast<char*>(&physical_time), sizeof(real));
    is.close();

    return physical_time;
  }


  void load_scalars(real* data, const std::string output_directory, const int coords[]) {

    std::string filepath = output_directory + "/restart/restart" +
      "_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) +
      "_" + std::to_string(coords[2]) + ".bin";

    int nscalars = 0;

    std::ifstream is(filepath.c_str(), std::ios::binary | std::ios::in);
    is.seekg(sizeof(int)+sizeof(real));
    is.read(reinterpret_cast<char*>(&nscalars), sizeof(int));
    is.read(reinterpret_cast<char*>(data), nscalars*sizeof(real));

    is.close();
  }


  void load_field(real* data, int field_number, const std::string output_directory, const int coords[]) {

    std::string filepath = output_directory + "/restart/restart" +
      "_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) +
      "_" + std::to_string(coords[2]) + ".bin";

    int nscalars = 0;
    int field_size = 0;

    std::ifstream is(filepath.c_str(), std::ios::binary | std::ios::in);
    is.seekg(sizeof(int)+sizeof(real));
    is.read(reinterpret_cast<char*>(&nscalars), sizeof(int));
    is.seekg(nscalars*sizeof(real), std::ios::cur);

    for (int i = 0; i <= field_number; ++i) {

      is.read(reinterpret_cast<char*>(&field_size), sizeof(int));

      if (i == field_number) {
        is.read(reinterpret_cast<char*>(data), field_size*sizeof(real));
        break;
      } else {
        is.seekg(field_size*sizeof(real), std::ios::cur);
      }
    }

    is.close();
  }

} // namespace restart


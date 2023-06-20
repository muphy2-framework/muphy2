/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "utilities.h"
#include <string>
#include <unordered_map>
#include <iostream>
#include <cstdlib> // exit
#include "framework/parameter.h"

namespace util {

  void copy_array(real* target_array, const real* source_array, int array_length) {
    utilities_copy_array_(target_array, source_array, array_length);
  }

  void copy_array_no_gpu(real* target_array, const real* source_array, int array_length) {
    utilities_copy_array_no_gpu_(target_array, source_array, array_length);
  }

  void add_arrays(real* target_array, const real* array1, const real* array2,
                  real factor1, real factor2, int array_length) {
    utilities_add_arrays_(target_array, array1, array2, factor1, factor2, array_length);
  }

  void add_arrays(real* target_array, const real* array1, const real* array2, int array_length) {
    utilities_add_arrays_(target_array, array1, array2, 1., 1., array_length);
  }

  void array_add_equal(real* target_array, const real* source_array, int array_length) {
    utilities_add_arrays_(target_array, target_array, source_array, 1., 1., array_length);
  }

  void array_add_scalar(real* array, const real& scalar, int array_length) {
    utilities_array_add_scalar_(array, scalar, array_length);
  }

  void array_fill(real* array, const real& fill_value, int array_length) {
    utilities_array_fill_(array, fill_value, array_length);
  }

  void vector_magnitude(real* magnitude, const real* vector, const Parameter& p) {
    utilities_vector_magnitude_(magnitude, vector, p.res_x_minus_one, p.bd);
  }

  void str_cpp_to_fortran(const std::string cpp_string, char* fortran_string) {
    for (int i = 0; i < (int) cpp_string.length(); ++i)
      fortran_string[i] = cpp_string[i];
  }

  void check_umap_keys(const std::string keys[], const int nkeys, const Parameter& p) {
    
    // check for duplicate keys
    for (int i = 0; i < nkeys-1; ++i) {
      for (int j = i+1; j < nkeys; ++j) {
        if (keys[i] == keys[j]) {
          std::cerr<<"Error: Duplicate keys in parameter::setup_var/parameter::setup_var_bool/parameter::setup_var_int."<<std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
    }

    // check if key is present
    for (int i = 0; i < nkeys; ++i) {
      if (p.setup_var.find(keys[i]) == p.setup_var.end() &&
          p.setup_var_bool.find(keys[i]) == p.setup_var_bool.end() &&
          p.setup_var_int.find(keys[i]) == p.setup_var_int.end()) {

          std::cerr<<"Error: Key \"" << keys[i] <<
            "\" not present in any of parameter::setup_var, parameter::setup_var_bool and parameter::setup_var_int."<<std::endl;
          std::exit(EXIT_FAILURE);
      }
    }
  }

}


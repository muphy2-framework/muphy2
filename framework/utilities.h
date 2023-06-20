/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <string>
#include "framework/definitions.h"

struct Parameter;

namespace util {

  // fortran functions
  extern "C" void utilities_copy_array_(real* target_array, const real* source_array, const int& array_length);
  extern "C" void utilities_copy_array_no_gpu_(real* target_array, const real* source_array, const int& array_length);
  extern "C" void utilities_add_arrays_(real* target_array, const real* array1, const real* array2,
                                        const real& factor1, const real&  factor2, const int&  array_length);
  extern "C" void utilities_array_add_scalar_(real* array, const real& scalar, const int& array_length);
  extern "C" void utilities_array_fill_(real* array, const real& fill_value, const int& array_length);
  extern "C" void utilities_vector_magnitude_(real* magnitude, const real* vector, const int* dimX, const int* BD);

  // wrappers for fortran functions
  void copy_array(real* target_array, const real* source_array, int array_length);
  void copy_array_no_gpu(real* target_array, const real* source_array, int array_length);
  void array_add_scalar(real* array, const real& scalar, int array_length);
  void array_fill(real* array, const real& fill_value, int array_length);
  void add_arrays(real* target_array, const real* array1, const real* array2,
                  real factor1, real factor2, int array_length);
  void add_arrays(real* target_array, const real* array1, const real* array2, int array_length);
  void array_add_equal(real* target_array, const real* source_array, int array_length);
  void vector_magnitude(real* magnitude, const real* vector, const Parameter& p);

  void str_cpp_to_fortran(const std::string cpp_string, char* fortran_string);
  void check_umap_keys(const std::string keys[], const int nkeys, const Parameter& p);
}

#endif // UTILITIES_H_

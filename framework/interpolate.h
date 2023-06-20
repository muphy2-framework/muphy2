/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include "framework/definitions.h"

struct Parameter;
struct Species;

namespace ipol {

  // fortran functions
  extern "C" void interpolate_yee_face_to_centered_(const real* u_face, real* u_centered, const int* dimX, const int* BD);
  extern "C" void interpolate_yee_edge_to_centered_(const real* u_edge, real* u_centered, const int* dimX, const int* BD);
  extern "C" void interpolate_yee_centered_to_face_(const real* u_centered, real* u_face, const int* dimX, const int* BD);
  extern "C" void interpolate_linear_extrapolation_to_boundary_f_(real* f, const int* dimX, const int* dimV,
                                                                  const int* BD, const int& direction);

  // wrappers for fortran functions
  void yee_face_to_centered(const real* u_face, real* u_centered, const Parameter& p);
  void yee_edge_to_centered(const real* u_edge, real* u_centered, const Parameter& p);
  void yee_centered_to_face(const real* u_centered, real* u_face, const Parameter& p);
  void linear_extrapolation_to_boundary_f(real* f, const Parameter& p, const Species& s, int direction);
}

#endif // INTERPOLATE_H_

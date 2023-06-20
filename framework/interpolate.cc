/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "interpolate.h"
#include "framework/definitions.h"
#include "framework/parameter.h"

namespace ipol {

  void yee_face_to_centered(const real* u_face, real* u_centered, const Parameter& p)
  {
    interpolate_yee_face_to_centered_(u_face, u_centered, p.res_x_minus_one, p.bd);
  }

  void yee_edge_to_centered(const real* u_edge, real* u_centered, const Parameter& p)
  {
    interpolate_yee_edge_to_centered_(u_edge, u_centered, p.res_x_minus_one, p.bd);
  }

  void yee_centered_to_face(const real* u_centered, real* u_face, const Parameter& p)
  {
    interpolate_yee_centered_to_face_(u_centered, u_face, p.res_x_minus_one, p.bd);
  }

  void linear_extrapolation_to_boundary_f(real* f, const Parameter& p, const Species& s, int direction)
  {
    interpolate_linear_extrapolation_to_boundary_f_(f, p.res_x_minus_one, s.res_v, p.bd, direction);
  }
}

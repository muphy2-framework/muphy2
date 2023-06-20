/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef CRITERION_POSITION_H_
#define CRITERION_POSITION_H_

#include "framework/criterion.h"

// criterion based on position within the simulation domain

struct Parameter;
struct Block_id;


namespace criterion {

class Criterion_position : public Criterion {
public:
  Criterion_position(int nused_schemes, int* scheme_hierarchy,
      int* base_point_x, int* base_point_y, int* base_point_z,
      int* width, int* height, int* depth);
  ~Criterion_position();
  int get_criterion_id();
  int evaluate(const Block_id& block_id, const Parameter& p);
private:
  int nused_schemes_;
  int* scheme_hierarchy_;

  // location of the area covered by a scheme (in MPI coordinates)
  // the area includes base_point and does not anymore include base_point+width
  int* base_point_x_,
     * base_point_y_,
     * base_point_z_,
     * width_,
     * height_,
     * depth_;
};

} // namespace criterion

#endif // CRITERION_POSITION_H_

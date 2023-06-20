/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "criterion_position.h"

#include "framework/parameter.h"
#include "framework/block.h"


namespace criterion {

Criterion_position::Criterion_position(int nused_schemes, int* scheme_hierarchy,
      int* base_point_x, int* base_point_y, int* base_point_z,
      int* width, int* height, int* depth) {

  nused_schemes_ = nused_schemes;

  scheme_hierarchy_ = new int[nused_schemes_];
  base_point_x_ = new int[nused_schemes_-1];
  base_point_y_ = new int[nused_schemes_-1];
  base_point_z_ = new int[nused_schemes_-1];
  width_ = new int[nused_schemes_-1];
  height_ = new int[nused_schemes_-1];
  depth_ = new int[nused_schemes_-1];

  for (int i = 0; i < nused_schemes; ++i) {
    scheme_hierarchy_[i] = scheme_hierarchy[i];
  }
  for (int i = 0; i < nused_schemes-1; ++i) {
    base_point_x_[i] = base_point_x[i];
    base_point_y_[i] = base_point_y[i];
    base_point_z_[i] = base_point_z[i];
    width_[i] = width[i];
    height_[i] = height[i];
    depth_[i] = depth[i];
  }
}


Criterion_position::~Criterion_position() {
  delete[] depth_;
  delete[] height_;
  delete[] width_;
  delete[] base_point_z_;
  delete[] base_point_y_;
  delete[] base_point_x_;
  delete[] scheme_hierarchy_;
}


int Criterion_position::get_criterion_id() { return Parameter::kCriterionPosition; }


int Criterion_position::evaluate(const Block_id& block_id, const Parameter& p) {

  for (int i = 0; i < nused_schemes_-1; ++i) {
    if (block_id.my_coords[0] >= base_point_x_[i] && block_id.my_coords[0] < base_point_x_[i] + width_[i] &&
        block_id.my_coords[1] >= base_point_y_[i] && block_id.my_coords[1] < base_point_y_[i] + height_[i] &&
        block_id.my_coords[2] >= base_point_z_[i] && block_id.my_coords[2] < base_point_z_[i] + depth_[i]) {
      return scheme_hierarchy_[i];
    }
  }
  return scheme_hierarchy_[nused_schemes_-1];

} // end evaluate()

} // namespace criterion


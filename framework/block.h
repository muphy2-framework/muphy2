/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef BLOCK_H_
#define BLOCK_H_

#include "framework/definitions.h"
#include "framework/mpi_boundary.h"
#include "framework/parameter.h"

namespace scheme{
  class Scheme;
}

struct Block_id {
  int my_id,
      neighbour_ids[6],
      my_coords[3],
      my_scheme_id,
      neighbour_scheme_ids[6],
      compute_node = -1,
      selected_device = -1;
};

class Block {

public:
  Block();
  ~Block();

  void loop_cycle(real &t);

  real get_t_begin();
  real get_t_end();

private:
  void update_scheme();
  void load_balancing();

  Block_id block_id_;
  Parameter parameter_;
  Mpi_boundary boundary_;
  scheme::Scheme *scheme_;

  MPI_Comm comm3d_;

  real dt_; // global time step size

  int output_number_ = 0,
      max_local_rank_;
  long int nsteps_total = 0;
  real dt_update_scheme_,
      dt_output_vtk_,
      dt_output_csv_,
      wall_time_total_s = 0.;
  bool freeze_dt_default_;
};

#endif // BLOCK_H_

/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cstdlib> // exit
//#include <cuda_runtime.h>
//#include <cuda_profiler_api.h>
#include "mpi.h"
#include "framework/definitions.h"
#include "framework/block.h"

int main(int argc, char** argv)
{
  MPI_Init(NULL,NULL);

  Block block;

  real t = block.get_t_begin(); // physical time

  //cudaProfilerStart();

  while(t < block.get_t_end())
  {
    block.loop_cycle(t);
  }

  //cudaDeviceSynchronize();
  //cudaProfilerStop();

  block.~Block();

  MPI_Finalize();

  std::exit(EXIT_SUCCESS);
}

/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef MHD_H_
#define MHD_H_

#include "framework/model.h"

// TODO this model is work in progress and not functional

struct Parameter;

// fortran functions
extern "C" void five_moment_mhd_single_step_rk3_(real* fluid, real* fluid_interim,
                                                 const real* E, const real* B,  const real* J,
                                                 const real& alpha0, const real& alpha1, const real& beta,
                                                 const real& cweno_limiter, const real& cweno_epsilon,
                                                 const int* dimX, const int* BD, const real* dx, const real& dt,
                                                 const int& dimensionality_x);

extern "C" real five_moment_mhd_max_dt_(const real* fluid, const real& cfl, const real& cweno_limiter,
                                        const int* dimX, const int* BD, const real* dx);

namespace model {

class MHD1Temperature: public Model {
public:
   MHD1Temperature(const Parameter& parameter);
  ~ MHD1Temperature();

  int get_model_id();
  real get_max_dt();

  // wrappers for fortran functions
  void single_step_rk3(const real dt, real* five_moments_tmp,
                       const real* E, const real* B, const real* J,
                       const real alpha0, const real alpha1, const real beta);

  real* five_moments;

private:
  const Parameter& parameter_;
};

} // namespace model

#endif // MHD_H_

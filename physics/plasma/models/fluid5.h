/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// -*- lsst-c++ -*-

/*
 * This file is part of the muphy2 framwork.
 *
 * Developed by members of Theoretical Physics I, Ruhr University Bochum
 *
 * TODO: LICENCE
 */

#ifndef FLUID5_H_
#define FLUID5_H_

#include "framework/model.h"

struct Parameter;
struct Species;

// fortran functions
extern "C" void five_moment_single_step_rk3_(real* fluid, real* fluid_interim, const real* E, const real* B,
                    const real& q_m, const real& alpha0, const real& alpha1, const real& beta, const real& cweno_limiter,
                    const real& cweno_epsilon, const int* dimX, const int* BD, const real* dx, const real& dt,
                    const int& dimensionality_x, const int& dimensionality_v);
extern "C" void five_moment_calc_un_from_flux_(real* un, const real* fluid, const real& cweno_limiter, const real& cweno_epsilon,
                    const int* dimX, const int* BD, const int& dimensionality_v);
extern "C" real five_moment_max_dt_(const real* fluid, const real& cfl, const real& cweno_limiter,
                    const int* dimX, const int* BD, const real* dx, const int& dimensionality_v);

namespace model {

class Fluid5 : public Model {
public:
  Fluid5(const Parameter& parameter, const Species& species);
  ~Fluid5();

  int get_model_id();
  real get_max_dt();

  // wrappers for fortran functions
  void single_step_rk3(const real dt, real* five_mom_tmp, const real* E, const real* B,
                       const real alpha0, const real alpha1, const real beta);
  void calc_un_from_flux(real* un);

  real* five_moments;

private:
  const Parameter& parameter_;
  const Species& species_;
};

} // namespace model

#endif // FLUID5_H_

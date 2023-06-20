/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef POISSON_H_
#define POISSON_H_

#include "framework/model.h"


struct Parameter;

// fortran functions
extern "C" void poisson_calc_e_from_phi_(real* E, const real* phi, const int* dimX,
                        const int* BD, const real* dx);
extern "C" void poisson_gauss_seidel_(real* phi, const real* rho, const real& eps0,
                        const int* dimX, const int* BD, const real* dx, const int& dimensionality_x);
extern "C" void poisson_sor_(real* phi, const real* rho, const real& eps0,
                        const int* dimX, const int* BD, const real* dx, const int& dimensionality_x);
extern "C" real poisson_difference_max_value_(const real* phi_old, const real* phi_new, const int* dimX, const int* BD);
extern "C" real poisson_sum_phi_(const real* phi, const int* dimX, const int* BD);

namespace model {

class Poisson : public Model {
public:
  Poisson(const Parameter& parameter);
  ~Poisson();

  int get_model_id();
  real get_max_dt();
  void solver_iterations(const real *rho);

  // wrappers for fortran functions
  void gauss_seidel(const real *rho);
  void sor(const real *rho);
  void calc_E_from_phi();
  real difference_max_value(const real *phi_old, const real *phi_new);
  real sum_phi();

  real* E;
  real* B;
  real* phi;

private:
  const Parameter &parameter_;

  real* phi_old_;
  const int iterations_per_convergence_check_ = 10;
};

} // namespace model

#endif // POISSON_H_

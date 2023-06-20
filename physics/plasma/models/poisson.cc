/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "poisson.h"
#include "framework/utilities.h"
#include "framework/parameter.h"

namespace model {

Poisson::Poisson(const Parameter &parameter) : parameter_(parameter) {
  E = new real[3 * parameter_.ncells_x]();
  B = new real[3 * parameter_.ncells_x]();
  phi = new real[parameter_.ncells_x]();
  phi_old_ = new real[parameter_.ncells_x]();

  #pragma acc enter data copyin(E[0:3*parameter_.ncells_x])
  #pragma acc enter data copyin(B[0:3*parameter_.ncells_x])
  #pragma acc enter data copyin(phi[0:parameter_.ncells_x])
  #pragma acc enter data copyin(phi_old_[0:parameter_.ncells_x])
}

Poisson::~Poisson() {
  #pragma acc exit data delete(phi_old_[0:parameter_.ncells_x])
  #pragma acc exit data delete(phi[0:parameter_.ncells_x])
  #pragma acc exit data delete(B[0:3*parameter_.ncells_x])
  #pragma acc exit data delete(E[0:3*parameter_.ncells_x])

  delete[] phi_old_;
  delete[] phi;
  delete[] B;
  delete[] E;
}


int Poisson::get_model_id() { return parameter_.kPoisson; }


real Poisson::get_max_dt() { return 1.e10; }


void Poisson::solver_iterations(const real *rho) {

  real max_difference = 100.*parameter_.poisson_convergence_threshold;

  while (max_difference > parameter_.poisson_convergence_threshold) {

    util::copy_array(phi_old_, phi, parameter_.ncells_x);
  
    for (int i = 0; i < iterations_per_convergence_check_; ++i) {
      if (parameter_.poisson_solver == "sor") {
        sor(rho);
      } else if (parameter_.poisson_solver == "gauss_seidel") {
        gauss_seidel(rho);
      }
    }

    max_difference = difference_max_value(phi_old_, phi);
  }
}


real Poisson::sum_phi() {
  const Parameter &p(parameter_);
  return poisson_sum_phi_(phi, p.res_x_minus_one, p.bd);
}


real Poisson::difference_max_value(const real *phi_old, const real *phi_new) {
  const Parameter &p(parameter_);
  return poisson_difference_max_value_(phi_old, phi_new, p.res_x_minus_one, p.bd);
}


void Poisson::gauss_seidel(const real *rho) {
  const Parameter &p(parameter_);
  poisson_gauss_seidel_(phi, rho, p.eps0, p.res_x_minus_one, p.bd, p.dx, p.dimensionality_x);
}


void Poisson::sor(const real *rho) {
  const Parameter &p(parameter_);
  poisson_sor_(phi, rho, p.eps0, p.res_x_minus_one, p.bd, p.dx, p.dimensionality_x);
}


void Poisson::calc_E_from_phi() {
  const Parameter &p(parameter_);
  poisson_calc_e_from_phi_(E, phi, p.res_x_minus_one, p.bd, p.dx);
}

} // namespace model

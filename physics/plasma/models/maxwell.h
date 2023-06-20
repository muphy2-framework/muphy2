/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef MAXWELL_H_
#define MAXWELL_H_

#include "framework/model.h"


struct Parameter;

// fortran functions
extern "C" void step_e_maxwell_fdtd_(real* E, const real* B, const real* j, const real& c, const real& mu0,
                        const int* dimX, const int* BD, const real* dx, const real& dt);
extern "C" void step_b_maxwell_fdtd_(const real* E, real* B,
                        const int* dimX, const int* BD, const real* dx, const real& dt);
extern "C" void step_b_ipol_maxwell_fdtd_(const real* E, real* B, real* B_ipol,
                        const int* dimX, const int* BD, const real* dx, const real& dt);

namespace model {

class Maxwell : public Model {
public:
  Maxwell(const Parameter& parameter);
  ~Maxwell();

  int get_model_id();
  real get_max_dt();

  // wrappers for fortran functions
  void step_E(const real dt, const real *j);
  void step_B(const real dt);
  void step_B_ipol(const real dt);

  real* E;
  real* B;
  real* B_ipol;   // magnetic field interpolated to integer times

  int substeps = 1;
  real dt = 0.;

private:
  const Parameter &parameter_;
};

} // namespace model

#endif // MAXWELL_H_

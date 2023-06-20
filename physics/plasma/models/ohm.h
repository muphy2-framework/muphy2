/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef OHM_H_
#define OHM_H_

#include "framework/model.h"

// TODO this model is work in progress and not functional

struct Parameter;

// fortran functions
extern "C" void ohm_calc_e_fc_ideal_(real* E_fc, const real* B_ec, const real* u_cc,
                                     const int* dimX, const int* BD);
extern "C" void ohm_calc_e_fc_hallterm_(real* E_fc, const real* B_ec, const real* J_fc,
				       const real* rho_cc, const int* dimX, const int* BD);
extern "C" void ohm_calc_e_fc_pressureterm_(real* E_fc, const real* P_cc, const real* rho_cc,
					    const real* dx, const int* dimX, const int* BD);
extern "C" void ohm_calc_e_fc_resistivityterm_(real* E_fc, const real* J_fc, const real& resistivity,
					       const int* dimX, const int* BD);
extern "C" void ohm_step_b_(const real* E, real* B,
                                  const int* dimX, const int* BD, const real* dx, const real& dt);
extern "C" void ohm_calc_j_(const real* J, real* B,
                                  const int* dimX, const int* BD, const real* dx);

namespace model {

class Ohm : public Model { //@TODO: rename OhmBase ?
 public:
  Ohm(const Parameter& parameter);
  ~Ohm();

  int get_model_id();
  real get_max_dt();

  // wrappers for fortran functions
  void step_B(const real dt);
  void calc_J();

  real *E;
  real *B;
  real *J;

 protected:
  const Parameter &parameter_;
};

class OhmResistive : public Ohm {
 public:
  using Ohm::Ohm;
  void calc_E_fc(const real *u);
};

class OhmHall : public Ohm {
 public:
  using Ohm::Ohm;
  void calc_E_fc(const real *u, const real *rho);
};

class OhmPressure : public Ohm {
 public:
  using Ohm::Ohm;
  void calc_E_fc(const real *u, const real *P, const real *rho);
};

class OhmHallPressure : public Ohm {
 public:
  using Ohm::Ohm;
  void calc_E_fc(const real *u, const real *P, const real *rho);
};

} // namespace model

#endif // OHM_H_

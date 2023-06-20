/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef VEVIP_H_
#define VEVIP_H_

#include "framework/scheme.h"

struct Block_id;
class Mpi_boundary;
struct Parameter;

namespace model {
class Model;
class Vlasov;
class Poisson;
}

namespace scheme {
// Vlasov electrons, Vlasov ions, Poisson

class VeViP : public Scheme {
public:
  VeViP(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter);
  ~VeViP();

  int get_scheme_id();
  void init();
  void step();
  real get_dt(int species = 0);
  void set_dt(real dt) { dt_ = dt; };
  void output(real t, int output_number, bool output_vtk);
  int evaluate_criterion();

private:
  void exchange_plasma_solvers(real* data_e, real* data_i);

  const Block_id &block_id_;
  Mpi_boundary &boundary_;
  const Parameter &parameter_;

  model::Vlasov *vlasov_e_;
  model::Vlasov *vlasov_i_;
  model::Poisson *poisson_;

  real *rho_;

  real dt_ = 0.;
  const bool include_boundary_cells_ = true;
};

} // namespace scheme

#endif // VEVIP_H_

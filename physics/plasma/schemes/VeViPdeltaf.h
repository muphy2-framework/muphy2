/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef VEVIPDELTAF_H_
#define VEVIPDELTAF_H_

#include "framework/scheme.h"

struct Block_id;
class Mpi_boundary;
struct Parameter;

namespace model {
class Model;
class Vlasov;
class Fluid10;
class Poisson;
}

namespace scheme {
// Vlasov electrons, Vlasov ions, Poisson

class VeViPdeltaf : public Scheme {
public:
  VeViPdeltaf(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter);
  ~VeViPdeltaf();

  int get_scheme_id();
  void init();
  void step();
  real get_dt(int species = 0);
  void set_dt(real dt) { dt_ = dt; };
  void output(real t, int output_number, bool output_vtk);
  int evaluate_criterion();

private:
  void exchange_plasma_solvers(real* data_vlasov_e, real* data_vlasov_i,
          real* data_fluid_e, real* data_fluid_i, bool runge_kutta_only=false);

  const Block_id &block_id_;
  Mpi_boundary &boundary_;
  const Parameter &parameter_;

  model::Vlasov *vlasov_e_;
  model::Vlasov *vlasov_i_;
  model::Fluid10 *fluid_e_;
  model::Fluid10 *fluid_i_;
  model::Poisson *poisson_;

  real *rho_,
       *E_old_,
       *Q_raw_old_e_,
       *Q_raw_old_i_,
       *Q_raw_new_e_,
       *Q_raw_new_i_,
       *ten_moments_tmp_e_,
       *ten_moments_tmp_i_,
       *total_flux_e_,
       *total_flux_i_,
       *source_e_,
       *source_i_;

  real dt_ = 0.;
  const bool include_boundary_cells_ = true,
             runge_kutta_only_ = true;
};

} // namespace scheme

#endif // VEVIPDELTAF_H_

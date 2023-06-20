/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef F10EVIMDELTAF_H_
#define F10EVIMDELTAF_H_

#include "framework/scheme.h"

struct Block_id;
class Mpi_boundary;
struct Parameter;

namespace model {
class Model;
class Vlasov;
class Fluid10;
class Maxwell;
}

namespace scheme {

// ten moments fluid electrons, ten moments fluid ions, Maxwell

class F10eViMdeltaf : public Scheme {
public:
  F10eViMdeltaf(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter);
  ~F10eViMdeltaf();

  int get_scheme_id();
  void init();
  void step();
  real get_dt(int species = 0);
  void set_dt(real dt) { dt_ = dt; };
  void output(real t, int output_number, bool output_vtk);
  int evaluate_criterion();
  void replace_dealloc_old_models(model::Model* plasma_model[],
          model::Model* electromagnetic_model);
  void calc_alloc_converted_model(model::Model*& model_out,
          int species, int target_scheme_id);
  void send_model_data(real* out_buffer, int receiver_id);
  void receive_model_data(real* in_buffer, int sender_id);
  int get_data_size();

  // for electron subcycling
  int get_substeps() { return electron_substeps_; };
  int get_substep_counter() { return electron_substep_counter_; };
  void set_substeps(int substeps) { electron_substeps_ = substeps; };
  void set_substep_counter(int substep_counter) { electron_substep_counter_ = substep_counter; };

private:
  void heat_flux_closure(model::Fluid10 *fluid, const real *ten_moments_tmp, real *source);
  void exchange_plasma_solvers(real* data_vlasov_i,
          real* data_fluid_e, real* data_fluid_i, bool runge_kutta_only=false);

  const Block_id &block_id_;
  Mpi_boundary &boundary_;
  const Parameter &parameter_;

  model::Fluid10 *fluid_e_;
  model::Fluid10 *fluid_i_;
  model::Vlasov *vlasov_i_;
  model::Maxwell *maxwell_;

  real *un_e_,
       *un_i_,
       *j_fc_,        // Yee grid face centered
       *E_cc_,        // cell centered
       *B_cc_,
       *Q_raw_old_i_,
       *Q_raw_new_i_,
       *ten_moments_tmp_e_,
       *ten_moments_tmp_i_,
       *total_flux_e_,
       *total_flux_i_,
       *source_e_,
       *source_i_;

  real dt_ = 0.;
  const bool runge_kutta_only_ = true,
             include_boundary_cells_ = true;
  int electron_substeps_ = 1,
      electron_substep_counter_ = 0;
  bool ion_step_ = false;
};

} // namespace scheme

#endif // F10EVIMDELTAF_H_

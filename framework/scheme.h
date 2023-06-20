/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef SCHEME_H_
#define SCHEME_H_

#include "framework/definitions.h"

// scheme base class

namespace model {
class Model;
}

namespace scheme {

class Scheme {
public:
  virtual ~Scheme() { };
  virtual int get_scheme_id() = 0;
  virtual void init() = 0;
  virtual void step() = 0;
  virtual real get_dt(int species = 0) = 0;
  virtual void set_dt(real dt) = 0;
  virtual void output(real t, int output_number, bool vtk_output) = 0;
  virtual int evaluate_criterion() = 0;
  virtual void send_model_data(real* out_buffer, int receiver_id) { };
  virtual void receive_model_data(real* in_buffer, int sender_id) { };
  virtual int get_data_size() { return 0; }; // total number of model data array elements, needed for load balancing exchange
  virtual void replace_dealloc_old_models(model::Model* plasma_model[],
                  model::Model* electromagnetic_model) { };
  virtual void calc_alloc_converted_model(model::Model*& model_out,
                  int species, int target_scheme_id) { };
  virtual int get_substeps() { return 1; };
  virtual int get_substep_counter() { return 0; };
  virtual void set_substeps(int substeps) { };
  virtual void set_substep_counter(int substep_counter) { };
};

} // namespace scheme

#endif // SCHEME_H_

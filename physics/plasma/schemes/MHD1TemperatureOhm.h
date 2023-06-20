/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef MHDO_H_
#define MHDO_H_

#include "framework/scheme.h"

// TODO this scheme is work in progress and not functional

struct Block_id;
class Mpi_boundary;
struct Parameter;

namespace model {
class MHD1Temperature;
class OhmResistive;
class OhmHall;
class OhmPressure;
class OhmHallPressure;
}

namespace scheme {

// five moments fluid electrons, five moments fluid ions, Maxwell

class MHD1TemperatureOhm : public Scheme {
public:
  MHD1TemperatureOhm(const Block_id &block_id, Mpi_boundary &mpi_boundary, const Parameter &parameter);
  ~MHD1TemperatureOhm();

  int get_scheme_id();
  void init();
  void step();
  real get_dt(int species = 0);
  void set_dt(real dt) { dt_ = dt; };
  void output(real t, int output_number, bool output_vtk);
  int evaluate_criterion();

private:
  const Block_id &block_id_;
  Mpi_boundary &boundary_;
  const Parameter &parameter_;

  model::MHD1Temperature *mhd_;
  //model::OhmResistive *ohm_;
  //model::OhmHall *ohm_;
  model::OhmPressure *ohm_;
  //model::OhmHallPressure *ohm_;

  real *E_cc_,      // cell centered
       *B_cc_,
       *J_cc_,
       *five_moments_tmp_;

  real dt_ = 0.;
};

} // namespace scheme

#endif // MHDO_H_

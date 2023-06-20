/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "vlasov.h"
#include <iostream>
#include "framework/parameter.h"

namespace model {

Vlasov::Vlasov(const Parameter &parameter, const Species &species)
    : parameter_(parameter), species_(species) {
  f = new real[species_.ncells_xv]();
  #pragma acc enter data copyin(f[0:species_.ncells_xv])
}

Vlasov::~Vlasov() {
  #pragma acc exit data delete(f[0:species_.ncells_xv])
  delete[] f;
}

int Vlasov::get_model_id() { return parameter_.kVlasov; }

real Vlasov::get_max_dt() {
  const Parameter &p(parameter_);
  const Species &s(species_);

  real cfl = p.cfl;
  if (p.cfl_vlasov_ions > 0. && s.q > 0.) {
    cfl = p.cfl_vlasov_ions;
  }

  real dt = cfl * (p.xe[0] - p.xb[0]) / p.res_x_total[0] / s.ve[0];

  if (p.dimensionality_x > 1)
    dt = std::min(dt, cfl * (p.xe[1] - p.xb[1]) / p.res_x_total[1] / s.ve[1]);

  if (p.dimensionality_x > 2)
    dt = std::min(dt, cfl * (p.xe[2] - p.xb[2]) / p.res_x_total[2] / s.ve[2]);

  return dt;
}

void Vlasov::step_x(real dt) {
  const Parameter &p(parameter_);
  const Species &s(species_);

  if (p.vlasov_solver == "pfc") {
    step_x_vlasov_pfc_(f, p.res_x_minus_one, s.res_v, p.bd, p.dx, s.dv, s.vb, dt, p.dimensionality_x);
  } else if (p.vlasov_solver == "lagrange5") {
    step_x_vlasov_lagrange5_(f, p.res_x_minus_one, s.res_v, p.bd, p.dx, s.dv, s.vb, dt, p.dimensionality_x);
  } else {
    std::cerr<<std::endl<<"Error: Invalid value for parameter::vlasov_solver.";
  }
}

void Vlasov::step_y(real dt) {
  const Parameter &p(parameter_);
  const Species &s(species_);
  if (p.vlasov_solver == "pfc") {
    step_y_vlasov_pfc_(f, p.res_x_minus_one, s.res_v, p.bd, p.dx, s.dv, s.vb, dt, p.dimensionality_x);
  } else if (p.vlasov_solver == "lagrange5") {
    step_y_vlasov_lagrange5_(f, p.res_x_minus_one, s.res_v, p.bd, p.dx, s.dv, s.vb, dt, p.dimensionality_x);
  } else {
    std::cerr<<std::endl<<"Error: Invalid value for parameter::vlasov_solver.";
  }
}

void Vlasov::step_z(real dt) {
  const Parameter &p(parameter_);
  const Species &s(species_);
  if (p.vlasov_solver == "pfc") {
    step_z_vlasov_pfc_(f, p.res_x_minus_one, s.res_v, p.bd, p.dx, s.dv, s.vb, dt);
  } else if (p.vlasov_solver == "lagrange5") {
    step_z_vlasov_lagrange5_(f, p.res_x_minus_one, s.res_v, p.bd, p.dx, s.dv, s.vb, dt);
  } else {
    std::cerr<<std::endl<<"Error: Invalid value for parameter::vlasov_solver.";
  }
}

void Vlasov::step_vx(const real dt, const real &E, const real &B, bool include_boundary_cells) {
  const Parameter &p(parameter_);
  const Species &s(species_);
  if (p.vlasov_solver == "pfc") {
    step_vx_vlasov_pfc_(f, &E, &B, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
        s.q, s.m, dt, p.dimensionality_x, include_boundary_cells);
  } else if (p.vlasov_solver == "lagrange5") {
    step_vx_vlasov_lagrange5_(f, &E, &B, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
        s.q, s.m, dt, p.dimensionality_x, include_boundary_cells);
  } else {
    std::cerr<<std::endl<<"Error: Invalid value for parameter::vlasov_solver.";
  }
}

void Vlasov::step_vy(const real dt, const real &E, const real &B, bool include_boundary_cells) {
  const Parameter &p(parameter_);
  const Species &s(species_);
  if (p.vlasov_solver == "pfc") {
    step_vy_vlasov_pfc_(f, &E, &B, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
        s.q, s.m, dt, p.dimensionality_x, include_boundary_cells);
  } else if (p.vlasov_solver == "lagrange5") {
    step_vy_vlasov_lagrange5_(f, &E, &B, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
        s.q, s.m, dt, p.dimensionality_x, include_boundary_cells);
  } else {
    std::cerr<<std::endl<<"Error: Invalid value for parameter::vlasov_solver.";
  }
}

void Vlasov::step_vz(const real dt, const real &E, const real &B, bool include_boundary_cells) {
  const Parameter &p(parameter_);
  const Species &s(species_);
  if (p.vlasov_solver == "pfc") {
    step_vz_vlasov_pfc_(f, &E, &B, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
        s.q, s.m, dt, p.dimensionality_x, include_boundary_cells);
  } else if (p.vlasov_solver == "lagrange5") {
    step_vz_vlasov_lagrange5_(f, &E, &B, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
        s.q, s.m, dt, p.dimensionality_x, include_boundary_cells);
  } else {
    std::cerr<<std::endl<<"Error: Invalid value for parameter::vlasov_solver.";
  }
}

} // namespace model

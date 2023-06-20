/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef VLASOV_H_
#define VLASOV_H_

#include "framework/model.h"


struct Parameter;
struct Species;

// fortran functions
// positive flux conservative solver
extern "C" void step_x_vlasov_pfc_(real* f, const int* dimX, const int* dimV, const int* BD,
                        const real* dx, const real* dv, const real* v_b, const real& dt, const int& dimensionality_x);
extern "C" void step_y_vlasov_pfc_(real* f, const int* dimX, const int* dimV, const int* BD,
                        const real* dx, const real* dv, const real* v_b, const real& dt, const int& dimensionality_x);
extern "C" void step_z_vlasov_pfc_(real* f, const int* dimX, const int* dimV, const int* BD,
                        const real* dx, const real* dv, const real* v_b, const real& dt);
extern "C" void step_vx_vlasov_pfc_(real* f, const real* E, const real* B, const int* dimX,
                        const int* dimV, const int* BD, const real* dv, const real* v_b,
                        const real& q, const real& m, const real& dt, const int& dimensionality_x, const int& include_boundary_cells);
extern "C" void step_vy_vlasov_pfc_(real* f, const real* E, const real* B, const int* dimX,
                        const int* dimV, const int* BD, const real* dv, const real* v_b,
                        const real& q, const real& m, const real& dt, const int& dimensionality_x, const int& include_boundary_cells);
extern "C" void step_vz_vlasov_pfc_(real* f, const real* E, const real* B, const int* dimX,
                        const int* dimV, const int* BD, const real* dv, const real* v_b,
                        const real& q, const real& m, const real& dt, const int& dimensionality_x, const int& include_boundary_cells);
// lagrange polynomial solver
extern "C" void step_x_vlasov_lagrange5_(real* f, const int* dimX, const int* dimV, const int* BD,
                        const real* dx, const real* dv, const real* v_b, const real& dt, const int& dimensionality_x);
extern "C" void step_y_vlasov_lagrange5_(real* f, const int* dimX, const int* dimV, const int* BD,
                        const real* dx, const real* dv, const real* v_b, const real& dt, const int& dimensionality_x);
extern "C" void step_z_vlasov_lagrange5_(real* f, const int* dimX, const int* dimV, const int* BD,
                        const real* dx, const real* dv, const real* v_b, const real& dt);
extern "C" void step_vx_vlasov_lagrange5_(real* f, const real* E, const real* B, const int* dimX,
                        const int* dimV, const int* BD, const real* dv, const real* v_b,
                        const real& q, const real& m, const real& dt, const int& dimensionality_x, const int& include_boundary_cells);
extern "C" void step_vy_vlasov_lagrange5_(real* f, const real* E, const real* B, const int* dimX,
                        const int* dimV, const int* BD, const real* dv, const real* v_b,
                        const real& q, const real& m, const real& dt, const int& dimensionality_x, const int& include_boundary_cells);
extern "C" void step_vz_vlasov_lagrange5_(real* f, const real* E, const real* B, const int* dimX,
                        const int* dimV, const int* BD, const real* dv, const real* v_b,
                        const real& q, const real& m, const real& dt, const int& dimensionality_x, const int& include_boundary_cells);

namespace model {

class Vlasov : public Model {
public:
  Vlasov(const Parameter& parameter, const Species& species);
  ~Vlasov();

  int get_model_id();
  real get_max_dt();

  // wrappers for fortran functions
  void step_x(const real dt);
  void step_y(const real dt);
  void step_z(const real dt);
  void step_vx(const real dt, const real& E, const real& B, bool include_boundary_cells=false);
  void step_vy(const real dt, const real& E, const real& B, bool include_boundary_cells=false);
  void step_vz(const real dt, const real& E, const real& B, bool include_boundary_cells=false);

  real* f;

private:
  const Parameter& parameter_;
  const Species& species_;
};

} // namespace model

#endif // VLASOV_H_

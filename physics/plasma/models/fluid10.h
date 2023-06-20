/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef FLUID10_H_
#define FLUID10_H_

#include "framework/model.h"

/*
    fluid(:,:,:,1) = n
    fluid(:,:,:,2:4) = n u
    fluid(:,:,:,5:10) = Eps_ij = \int v_i v_j f dv (energy density tensor)
    eps = 3 n k T / m + n u \dot u (scalar energy density)
    P_ij = m \int (v_i - u_i) (v_j - u_j) f dv (pressure tensor)
         = m Eps_ij - m n u_i u_j
    p = tr(p_ij)/N (scalar pressure; N dimensionality)
*/

struct Parameter;
struct Species;

// fortran functions
extern "C" void ten_moment_calc_single_step_rk3_(real* fluid_interim, real* total_flux, real* source,
                    const real* E, const real* B, const real& q_m, const real& cweno_limiter, const real& cweno_epsilon,
                    const int& dimensionality_x, const int* dimX, const int* BD, const real* dx);

extern "C" void ten_moment_deltaf_calc_single_step_rk3_(real* fluid_interim, real* total_flux, real* source, const real* Q_raw,
                    const real* E, const real* B, const real& q_m, const real& cweno_limiter, const real& cweno_epsilon,
                    const int& dimensionality_x, const int* dimX, const int* BD, const real* dx);

extern "C" void ten_moment_apply_single_step_rk3_(real* fluid, real* fluid_interim, real* total_flux, real* source,
                    const real& alpha0, const real& alpha1, const real& beta, const real& cweno_limiter,
                    const int& dimensionality_v, const int* dimX, const int* BD, const real& dt);

extern "C" void ten_moment_add_divheatflux_to_source_(const real* divQ_h, real* source, const real& m, const int* dimX);

extern "C" void ten_moment_isotropization_heat_flux_closure_(const real* fluid_interim, real* source, const real& k0,
                    const real& cweno_limiter, const int& dimensionality_v, const int* dimX, const int* BD);

extern "C" void ten_moment_gradient_p_heat_flux_closure_cycle_(const real* fluid_interim, real* source, real* P_m_interim,
                    const real& k0, const int& current_cycle, const int& ncycles, const real& cweno_limiter,
                    const int* dimX, const int* BD, const real* dx, const int& dimensionality_x, const int& dimensionality_v,
                    const real& dt);

extern "C" void ten_moment_gradient_t_heat_flux_closure_cycle_(const real* fluid_interim, real* source, real* T_m_interim,
                    const real& k0, const int& current_cycle, const int& ncycles, const real& cweno_limiter,
                    const int* dimX, const int* BD, const real* dx, const int& dimensionality_x, const real& dt);

extern "C" void ten_moment_gradient_sym_heat_flux_closure_cycle_(const real* fluid_interim, real* source, real* T_m_interim,
                    const real& k0, const int& current_cycle, const int& ncycles, const real& cweno_limiter,
                    const int* dimX, const int* BD, const real* dx, const real& dt);

extern "C" void ten_moment_calc_un_from_flux_(real* un, const real* fluid, const real& cweno_limiter, const real& cweno_epsilon,
                    const int* dimX, const int* BD);

extern "C" real ten_moment_max_dt_(const real* fluid, const real& cfl, const real& cweno_limiter,
                    const int* dimX, const int* BD, const real* dx);

namespace model {

class Fluid10 : public Model {
public:
  Fluid10(const Parameter& parameter, const Species& species);
  ~Fluid10();

  int get_model_id();
  real get_max_dt();

  // wrappers for fortran functions
  void calc_single_step_rk3(real* ten_moments_tmp, real* total_flux, real* source, const real* E, const real* B);
  void calc_single_step_rk3(real* ten_moments_tmp, real* total_flux, real* source, const real* Q_raw, const real* E, const real* B);
  void apply_single_step_rk3(const real dt, real* ten_moments_tmp, real* total_flux, real* source,
                             const real alpha0, const real alpha1, const real beta);
  void add_divheatflux_to_source(const real* divQ_h, real* source);
  void isotropization_heat_flux_closure(const real* ten_moments_tmp, real* source);
  void gradient_heat_flux_closure_cycle(const int iteration, const real dt,
                                        const real* ten_moments_tmp, real* source, real* P_m_tmp);
  void calc_un_from_flux(real* un, bool deltaf_scheme=false);

  real* ten_moments;

private:
  const Parameter& parameter_;
  const Species& species_;
};

} // namespace model

#endif // FLUID10_H_

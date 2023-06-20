/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "converters.h"
#include "framework/definitions.h"
#include "framework/parameter.h"
#include "framework/block.h"

namespace convert {

//_______________ converters for physical quantities _______________
void f_to_n(real *n, const real *f, const Species &s, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = true;
  convert_f_to_n_(n, f, array_size, s.res_v, s.dv, p.dimensionality_v, p.bd, skip_boundary);
}

void f_to_rho(real *rho, const real *f_e, const real *f_i, const Parameter &p) {
  convert_f_to_rho_(rho, f_e, f_i, p.species[p.kElectron].q, p.species[p.kIon].q, p.res_x_minus_one,
                  p.species[p.kElectron].res_v, p.species[p.kIon].res_v, p.bd,
                  p.species[p.kElectron].dv, p.species[p.kIon].dv, p.dimensionality_v);
}

void f_to_u(real *u, const real *f, const real *n, const Species &s, const Parameter &p) {
  convert_f_to_u_(u, f, n, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v);
}

void f_to_j(real *j, const real *f_e, const real *f_i, const Parameter &p) {
  convert_f_to_j_(j, f_e, f_i, p.species[p.kElectron].q, p.species[p.kIon].q, p.res_x_minus_one,
                  p.species[p.kElectron].res_v, p.species[p.kIon].res_v, p.bd, p.species[p.kElectron].dv,
                  p.species[p.kIon].dv, p.species[p.kElectron].vb, p.species[p.kIon].vb, p.dimensionality_v);
}

void f_to_un(real *un, const real *f, const Species &s, const Parameter &p) {
  convert_f_to_un_(un, f, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v);
}

void f_to_eps(real *eps, const real *f, const Species &s, const Parameter &p) {
  convert_f_to_eps_(eps, f, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v);
}

void f_to_heatflux(real *Q_h, const real *f, const real *u, const Species &s, const Parameter &p) {
  convert_f_to_heatflux_(Q_h, f, u, s.m, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v);
}

void f_to_raw_heatflux(real *Q_raw, const real *f, const Species &s,
                       const bool *exchange_in_this_direction, const Parameter &p) {
  int exchange_in_this_direction_f[6];
  for (int i = 0; i < 6; ++i) {
    exchange_in_this_direction_f[i] = (int) exchange_in_this_direction[i];
  }
  convert_f_to_raw_heatflux_(Q_raw, f, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb,
                             p.dimensionality_v, exchange_in_this_direction_f);
}

void u_to_j(real *j, const real *u_e, const real *u_i, const real *n_e,
            const real *n_i, const Parameter &p) {
  convert_u_to_j_(j, u_e, u_i, n_e, n_i, p.res_x_minus_one, p.bd, p.species[p.kElectron].q,
                  p.species[p.kIon].q);
}

void un_to_j(real *j, const real *un_e, const real *un_i, const Parameter &p) {
  convert_un_to_j_(j, un_e, un_i, p.res_x_minus_one, p.bd, p.species[p.kElectron].q,
                   p.species[p.kIon].q);
}

void un_to_j(real *j_face, const real *un_e_face, const real *un_i_face,
             const real *un_e_center, const real *un_i_center,
             const Parameter &p, const Block_id &id) {
  // special treatment for reflecting wall boundary conditions
  
  int reflecting_wall[6] = {0, 0, 0, 0, 0, 0};
  
  if (id.my_coords[0] == 0 && p.bd_cond_v[0][p.kXL] == 'z') {
    reflecting_wall[p.kXL] = 1;
  }
  if (id.my_coords[0] == p.nproc[0]-1 && p.bd_cond_v[0][p.kXU] == 'z') {
    reflecting_wall[p.kXU] = 1;
  }
  if (id.my_coords[1] == 0 && p.bd_cond_v[1][p.kYL] == 'z') {
    reflecting_wall[p.kYL] = 1;
  }
  if (id.my_coords[1] == p.nproc[1]-1 && p.bd_cond_v[1][p.kYU] == 'z') {
    reflecting_wall[p.kYU] = 1;
  }
  if (id.my_coords[2] == 0 && p.bd_cond_v[2][p.kZL] == 'z') {
    reflecting_wall[p.kZL] = 1;
  }
  if (id.my_coords[2] == p.nproc[2]-1 && p.bd_cond_v[2][p.kZU] == 'z') {
    reflecting_wall[p.kZU] = 1;
  }

  convert_un_to_j_reflecting_wall_(j_face, un_e_face, un_i_face, un_e_center, un_i_center,
                   p.res_x_minus_one, p.bd, p.species[p.kElectron].q, p.species[p.kIon].q,
                   reflecting_wall);
}

void un_to_u(real *u, const real *un, const real *n, const Parameter &p) {
  convert_un_to_u_(u, un, n, p.res_x_minus_one, p.bd);
}

void heatflux_to_divheatflux(real *divQ_h, const real *Q_h, const Parameter &p) {
  convert_heatflux_to_divheatflux_(divQ_h, Q_h, p.dx, p.res_x_minus_one, p.bd);
}

void mhd_1temperature_to_electron_pressure(real *P_e, const real *five_moments, const Parameter &p) {
  convert_mhd_1temperature_to_electron_pressure_(P_e, five_moments, p.species[p.kElectron].m, p.species[p.kElectron].T0,
                                                 p.species[p.kIon].T0, p.res_x_minus_one, p.bd);
}


//_______________ setup converters _______________
void fluid_quantities_to_f(real *f, const real *n, const real *u, const real *T,
                           const Species &s, const Parameter &p) {
  convert_fluid_quantities_to_f_(f, n, u, T, s.m, p.res_x_minus_one, s.res_v,
                                 p.bd, s.dv, s.vb, p.dimensionality_v);
}

void fluid_quantities_to_f_separate_background(real *f, const real *n, const real *u, const real *T,
          const real *n_bg, const real *u_bg, const real *T_bg, const Species &s, const Parameter &p) {
  convert_fluid_quantities_to_f_separate_background_(f, n, u, T, n_bg, u_bg, T_bg, s.m, p.res_x_minus_one, s.res_v,
                                 p.bd, s.dv, s.vb, p.dimensionality_v);
}

void fluid_quantities_to_mhd(real *five_moments, const real *n_e, const real *n_i,
                             const real *u_e, const real *u_i, const real *T_e, const real *T_i,
                             const real m_e, const real m_i, const Parameter &p) {
  convert_fluid_quantities_to_mhd_(five_moments, n_e, n_i, u_e, u_i, T_e, T_i, m_e, m_i,
                                   p.res_x_minus_one, p.bd);
}

void fluid_quantities_to_five_moments(real *five_moments, const real *n,
                                      const real *u, const real *T, const real m,
                                      const Parameter &p) {
  convert_fluid_quantities_to_five_moments_(five_moments, n, u, T, m,
                                            p.res_x_minus_one, p.bd, p.dimensionality_v);
}

void fluid_quantities_to_ten_moments(real *ten_moments, const real *n,
                                     const real *u, const real *T, const real m,
                                     const Parameter &p) {
  convert_fluid_quantities_to_ten_moments_(ten_moments, n, u, T, m,
                                           p.res_x_minus_one, p.bd, p.dimensionality_v);
}

//_______________ output converters _______________
void f_to_output(const real *f, real *n, real *u, real *T, real *P,
                 const Species &s, const Parameter &p) {
  convert_f_to_output_(f, n, u, T, P, s.m, p.res_x_minus_one, s.res_v, p.bd,
                       s.dv, s.vb, p.dimensionality_v);
}

void five_moments_to_output(const real *five_moments, real *n, real *u, real *T, real *P,
                            const real m, const Parameter &p) {
  convert_five_moments_to_output_(five_moments, n, u, T, P, m, p.res_x_minus_one,
                                  p.bd, p.dimensionality_v);
}

void ten_moments_to_output(const real *ten_moments, real *n, real *u, real *T,
                           real *P, const real m, const Parameter &p) {
  convert_ten_moments_to_output_(ten_moments, n, u, T, P, m, p.res_x_minus_one,
                                 p.bd, p.dimensionality_v);
}

//_______________ model converters _______________
void five_moments_to_ten_moments(real *ten_moments, const real *five_moments, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = false;
  convert_five_moments_to_ten_moments_(ten_moments, five_moments, array_size, p.dimensionality_v, p.bd, skip_boundary);
}

void five_moments_to_ten_moments(real *ten_moments, const real *five_moments, const int *array_size, const Parameter &p) {
  bool skip_boundary = false;
  convert_five_moments_to_ten_moments_(ten_moments, five_moments, array_size, p.dimensionality_v, p.bd, skip_boundary);
}

void ten_moments_to_five_moments(real *five_moments, const real *ten_moments, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = false;
  convert_ten_moments_to_five_moments_(five_moments, ten_moments, array_size, p.bd, skip_boundary);
}

void ten_moments_to_five_moments(real *five_moments, const real *ten_moments, const int *array_size, const Parameter &p) {
  bool skip_boundary = false;
  convert_ten_moments_to_five_moments_(five_moments, ten_moments, array_size, p.bd, skip_boundary);
}

void five_moments_to_f(real *f, const real *five_moments, const Species &s, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = false;
  convert_five_moments_to_f_(f, five_moments, s.m, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void five_moments_to_f(real *f, const real *five_moments, const int *array_size, const Species &s, const Parameter &p) {
  bool skip_boundary = false;
  convert_five_moments_to_f_(f, five_moments, s.m, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void ten_moments_to_f(real *f, const real *ten_moments, const Species &s, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = false;
  convert_ten_moments_to_f_(f, ten_moments, s.m, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void ten_moments_to_f(real *f, const real *ten_moments, const int *array_size, const Species &s, const Parameter &p) {
  bool skip_boundary = false;
  convert_ten_moments_to_f_(f, ten_moments, s.m, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void f_to_five_moments(real *five_moments, const real *f, const Species &s, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = false;
  convert_f_to_five_moments_(five_moments, f, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void f_to_five_moments(real *five_moments, const real *f, const int *array_size, const Species &s, const Parameter &p) {
  bool skip_boundary = false;
  convert_f_to_five_moments_(five_moments, f, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void f_to_ten_moments(real *ten_moments, const real *f, const Species &s, const Parameter &p) {
  int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
  bool skip_boundary = false;
  convert_f_to_ten_moments_(ten_moments, f, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void f_to_ten_moments(real *ten_moments, const real *f, const int *array_size, const Species &s, const Parameter &p) {
  bool skip_boundary = false;
  convert_f_to_ten_moments_(ten_moments, f, array_size, s.res_v, s.dv, s.vb, p.dimensionality_v, p.bd, skip_boundary);
}

void delta_f_fit(real *f, const real *n, const int *array_size, const Species &s, const Parameter &p, int direction) {
  if ((p.dimensionality_x == 2 && direction >= 4) || (p.dimensionality_x == 1 && direction >= 2)) {
    return;
  }
  convert_delta_f_fit_(f, n, array_size, s.m, p.res_x_minus_one, s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v, direction);
}

void delta_f_correction(real *f, real *ten_moments, const Species &s, const Parameter &p) {
  if (!p.cweno_calc_j_from_flux) {
    convert_delta_f_correction_keep_vlasov_n_(f, ten_moments, s.m, p.res_x_minus_one,
                                      s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v);
  } else {
    convert_delta_f_correction_(f, ten_moments, s.m, p.res_x_minus_one,
                                      s.res_v, p.bd, s.dv, s.vb, p.dimensionality_v);
  }
}

void five_ten_moments_fit(real *ten_moments, const Parameter &p, int direction) {
  if (p.five_ten_moment_coupling_fit) {
    convert_five_ten_moments_fit_(ten_moments, p.res_x_minus_one, p.bd, direction);
  }
}

} // namespace convert

/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef CONVERTERS_H_
#define CONVERTERS_H_

#include "framework/definitions.h"

struct Parameter;
struct Block_id;
struct Species;

namespace convert {

// fortran functions
//_______________ converters for physical quantities _______________
extern "C" void convert_f_to_n_(real *n, const real *f, const int *array_size,
                                const int *dimV, const real *dv, const int& dimensionality_v,
                                const int *BD, const int &skip_boundary);
extern "C" void convert_f_to_rho_(real *rho, const real *f_e, const real *f_i, const real &q_e, const real &q_i,
                                  const int *dimX, const int *dimV_e, const int *dimV_i, const int *BD,
                                  const real *dv_e, const real *dv_i, const int& dimensionality_v);
extern "C" void convert_f_to_u_(real *u, const real *f, const real *n, const int *dimX, const int *dimV,
                                const int *BD, const real *dv, const real *vb, const int& dimensionality_v);
extern "C" void convert_f_to_j_(real *j, const real *f_e, const real *f_i, const real &q_e, const real &q_i,
                                const int *dimX, const int *dimV_e, const int *dimV_i, const int *BD, const real *dv_e, const real *dv_i,
                                const real *vb_e, const real *vb_i, const int& dimensionality_v);
extern "C" void convert_f_to_un_(real *un, const real *f, const int *dimX, const int *dimV, const int *BD,
                                 const real *dv, const real *vb, const int& dimensionality_v);
extern "C" void convert_f_to_eps_(real *eps, const real *f, const int *dimX, const int *dimV, const int *BD,
                                  const real *dv, const real *vb, const int& dimensionality_v);
extern "C" void convert_f_to_heatflux_(real *Q_h, const real *f, const real *u, const real &m, const int *dimX, const int *dimV,
                                       const int *BD, const real *dv, const real *vb, const int &dimensionality_v);
extern "C" void convert_f_to_raw_heatflux_(real *Q_raw, const real *f, const int *dimX, const int *dimV, const int *BD, const real *dv,
                                           const real *vb, const int &dimensionality_v,
                                           const int *exchange_in_this_direction);
extern "C" void convert_u_to_j_(real *j, const real *u_e, const real *u_i, const real *n_e, const real *n_i,
                                const int *dimX, const int *BD, const real &q_e, const real &q_i);
extern "C" void convert_un_to_j_(real *j, const real *un_e, const real *un_i, const int *dimX, const int *BD, const real &q_e, const real &q_i);
extern "C" void convert_un_to_j_reflecting_wall_(real *j_face, const real *un_e_face, const real *un_i_face,
                                const real *un_e_center, const real *un_i_center, const int *dimX, const int *BD,
                                const real &q_e, const real &q_i, const int *reflecting_wall);
extern "C" void convert_un_to_u_(real *u, const real *un, const real *n, const int *dimX, const int *BD);
extern "C" void convert_heatflux_to_divheatflux_(real *divQ_h, const real *Q_h, const real *dx, const int *dimX, const int *BD);
extern "C" void convert_mhd_1temperature_to_electron_pressure_(real *P_e, const real *five_moments, const real &m_e,
                                                               const real &T_e, const real &T_i, const int *dimX, const int *BD);
//_______________ setup converters _______________
extern "C" void convert_fluid_quantities_to_f_(real *f, const real *n, const real *u, const real *T, const real &m,
                                               const int *dimX, const int *dimV, const int *BD, const real *dv,
                                               const real *vb, const int& dimensionality_v);
extern "C" void convert_fluid_quantities_to_f_separate_background_(real *f, const real *n, const real *u, const real *T,
                                               const real *n_bg, const real *u_bg, const real *T_bg, const real &m,
                                               const int *dimX, const int *dimV, const int *BD, const real *dv,
                                               const real *vb, const int& dimensionality_v);
extern "C" void convert_fluid_quantities_to_mhd_(real *five_moments, const real *n_e, const real *n_i,
						                                    const real *u_e, const real *u_i, const real *T_e, const real *T_i,
						                                    const real &m_e, const real &m_i, const int *dimX, const int *BD);
extern "C" void convert_fluid_quantities_to_five_moments_(real *five_moments, const real *n, const real *u, const real *T,
                                                          const real &m, const int *dimX, const int *BD, const int &dimensionality_v);
extern "C" void convert_fluid_quantities_to_ten_moments_(real *ten_moments, const real *n, const real *u, const real *T,
                                                         const real &m, const int *dimX, const int *BD, const int &dimensionality_v);
//_______________ output converters _______________
extern "C" void convert_f_to_output_(const real *f, real *n, real *u, real *T, real *P, const real &m,
                                     const int *dimX, const int *dimV, const int *BD,
                                     const real *dv, const real *vb, const int& dimensionality_v);
extern "C" void convert_five_moments_to_output_(const real *five_moments, real *n, real *u, real *T, real *P,
                                                const real &m, const int *dimX, const int *BD, const int& dimensionality_v);
extern "C" void convert_ten_moments_to_output_(const real *ten_moments, real *n, real *u, real *T, real *P,
                                               const real &m, const int *dimX, const int *BD, const int& dimensionality_v);
//_______________ model converters _______________
extern "C" void convert_ten_moments_to_five_moments_(real *five_moments, const real *ten_moments, const int *array_size,
                                                     const int *BD, const int &skip_boundary);
extern "C" void convert_five_moments_to_ten_moments_(real *ten_moments, const real *five_moments, const int *array_size,
                                                     const int& dimensionality_v, const int *BD, const int &skip_boundary);
extern "C" void convert_five_moments_to_f_(real *f, const real *five_moments, const real &m, const int *array_size,
                                          const int *dimV, const real *dv, const real *vb, const int& dimensionality_v,
                                          const int *BD, const int &skip_boundary);
extern "C" void convert_ten_moments_to_f_(real *f, const real *ten_moments, const real &m, const int *array_size,
                                         const int *dimV, const real *dv, const real *vb, const int& dimensionality_v,
                                         const int *BD, const int &skip_boundary);
extern "C" void convert_f_to_five_moments_(real *five_moments, const real *f, const int *array_size,
                                          const int *dimV, const real *dv, const real *vb, const int& dimensionality_v,
                                          const int *BD, const int &skip_boundary);
extern "C" void convert_f_to_ten_moments_(real *ten_moments, const real *f, const int *array_size,
                                         const int *dimV, const real *dv, const real *vb, const int& dimensionality_v,
                                         const int *BD, const int &skip_boundary);
extern "C" void convert_delta_f_fit_(real *f, const real *n, const int *array_size, const real &m, const int *dimX,
                                     const int *dimV, const int *BD, const real *dv, const real *vb,
                                     const int& dimensionality_v, const int& direction);
extern "C" void convert_delta_f_correction_(real *f, real *ten_moments, const real &m, const int *dimX,
                                     const int *dimV, const int *BD, const real *dv, const real *vb,
                                     const int& dimensionality_v);
extern "C" void convert_delta_f_correction_keep_vlasov_n_(real *f, real *ten_moments, const real &m, const int *dimX,
                                     const int *dimV, const int *BD, const real *dv, const real *vb,
                                     const int& dimensionality_v);
extern "C" void convert_five_ten_moments_fit_(real *ten_moments, const int *dimX, const int *BD, const int& direction);


// wrappers for fortran functions
//_______________ converters for physical quantities _______________
void f_to_n(real *n, const real *f, const Species &s, const Parameter &p);
void f_to_rho(real *rho, const real *f_e, const real *f_i, const Parameter &p);
void f_to_u(real *u, const real *f, const real *n, const Species &s, const Parameter &p);
void f_to_j(real *j, const real *f_e, const real *f_i, const Parameter &p);
void f_to_un(real *un, const real *f, const Species &s, const Parameter &p);
void f_to_eps(real *eps, const real *f, const Species &s, const Parameter &p);
void f_to_heatflux(real *Q_h, const real *f, const real *u, const Species &s, const Parameter &p);
void f_to_raw_heatflux(real *Q_raw, const real *f, const Species &s,
                       const bool *exchange_in_this_direction, const Parameter &p);
void u_to_j(real *j, const real *u_e, const real *u_i, const real *n_e, const real *n_i, const Parameter &p);
void un_to_j(real *j, const real *un_e, const real *un_i, const Parameter &p);
void un_to_j(real *j_face, const real *un_e_face, const real *un_i_face,
             const real *un_e_center, const real *un_i_center, const Parameter &p, const Block_id &id);
void un_to_u(real *u, const real *un, const real *n, const Parameter &p);
void heatflux_to_divheatflux(real *divQ_h, const real *Q_h, const Parameter &p);
void mhd_1temperature_to_electron_pressure(real *P_e, const real *five_moments, const Parameter &p);
//_______________ setup converters _______________
void fluid_quantities_to_f(real *f, const real *n, const real *u, const real *T, const Species &s, const Parameter &p);
void fluid_quantities_to_f_separate_background(real *f, const real *n, const real *u, const real *T,
                             const real *n_bg, const real *u_bg, const real *T_bg, const Species &s, const Parameter &p);
void fluid_quantities_to_mhd(real *five_moments, const real *n_e, const real *n_i,
                             const real *u_e, const real *u_i, const real *T_e, const real *T_i,
                             const real m_e, const real m_i, const Parameter &p);
void fluid_quantities_to_five_moments(real *five_moments, const real *n, const real *u, const real *T, const real m, const Parameter &p);
void fluid_quantities_to_ten_moments(real *ten_moments, const real *n, const real *u, const real *T, const real m, const Parameter &p);
//_______________ output converters _______________
void f_to_output(const real *f, real *n, real *u, real *T, real *P, const Species &s, const Parameter &p);
void five_moments_to_output(const real *five_moments, real *n, real *u, real *T, real *P, const real m, const Parameter &p);
void ten_moments_to_output(const real *ten_moments, real *n, real *u, real *T, real *P, const real m, const Parameter &p);
//_______________ model converters _______________
void five_moments_to_ten_moments(real *ten_moments, const real *five_moments, const Parameter &p);
void five_moments_to_ten_moments(real *ten_moments, const real *five_moments, const int *array_size, const Parameter &p);
void ten_moments_to_five_moments(real *five_moments, const real *ten_moments, const Parameter &p);
void ten_moments_to_five_moments(real *five_moments, const real *ten_moments, const int *array_size, const Parameter &p);
void five_moments_to_f(real *f, const real *five_moments, const Species &s, const Parameter &p);
void five_moments_to_f(real *f, const real *five_moments, const int *array_size, const Species &s, const Parameter &p);
void ten_moments_to_f(real *f, const real *ten_moments, const Species &s, const Parameter &p);
void ten_moments_to_f(real *f, const real *ten_moments, const int *array_size, const Species &s, const Parameter &p);
void f_to_five_moments(real *five_moments, const real *f, const Species &s, const Parameter &p);
void f_to_five_moments(real *five_moments, const real *f, const int *array_size, const Species &s, const Parameter &p);
void f_to_ten_moments(real *ten_moments, const real *f, const Species &s, const Parameter &p);
void f_to_ten_moments(real *ten_moments, const real *f, const int *array_size, const Species &s, const Parameter &p);
void delta_f_fit(real *f, const real *n, const int *array_size, const Species &s, const Parameter &p, int direction);
void delta_f_correction(real *f, real *ten_moments, const Species &s, const Parameter &p);
void five_ten_moments_fit(real *ten_moments, const Parameter &p, int direction);
}

#endif // CONVERTERS_H_

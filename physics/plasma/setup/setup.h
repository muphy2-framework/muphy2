/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef SETUP_H_
#define SETUP_H_

#include "framework/definitions.h"

struct Parameter;
namespace model {
class Model;
}

namespace setup {

// fortran functions
extern "C" void setup_landau_damping_(real *f_e, real *f_i, real *E, real *B, const real &m_e,
        const real &T0_e, const real &eps0, const real &alpha, const int *dimX, const int *dimV_e, const int *dimV_i,
        const int *BD, const real *dx, const real *vb_e, const real *dv_e, const real *dv_i,
        const int& dimensionality_x, const int& dimensionality_v, const real *xb_loc);

extern "C" void setup_two_stream_instability_(real *f_e, real *f_i, real *E, real *B, const real &m_e,
        const real &T0_e, const real &alpha, const int *dimX, const int *dimV_e, const int *dimV_i,
        const int *BD, const real *dx, const real *vb_e, const real *dv_e,
        const real *dv_i, const int& dimensionality_v, const real *xb_loc);

extern "C" void setup_orszag_tang_(real *n_e, real *n_i, real *u_e, real *u_i, real *E, real *B,
        const real &c, const real &L, const real &n_bg, const real &delta_u,
        const int *dimX, const int *BD, const real *dx, const real *xb_loc);

extern "C" void setup_orszag_tang_3d_(real *n_e, real *n_i, real *u_e, real *u_i, real *E, real *B,
        const real &L, const real &L_z, const real &n_bg, const real &delta_u,
        const real &cross_helicity, const int *dimX, const int *BD, const real *dx, const real *xb_loc);

extern "C" void setup_harris_sheet_(real *n_e, real *n_i, real *u_e, real *u_i, real *n_e_bg, real *n_i_bg,
        real *u_e_bg, real *u_i_bg, const real &T0_e, const real &T0_i, real *E, real *B, const real &lambda, const real &psi,
        const real &n_bg, const real &guide_field, const real &noise_level, const int &drifting_background, const int &sine_perturbation,
        const int *dimX, const int *BD, const real *xb, const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_double_harris_sheet_(real *n_e, real *n_i, real *u_e, real *u_i,
        const real &T0_e, const real &T0_i, real *E, real *B, const real &lambda, const real &psi,
        const real &n_bg, const real &guide_field, const real &noise_level, const int *dimX,
        const int *BD, const real *xb, const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_double_harris_sheet_gaussian_hump_(real *n_e, real *n_i, real *u_e, real *u_i, real *n_e_bg, real *n_i_bg,
        real *u_e_bg, real *u_i_bg, const real &T0_e, const real &T0_i, real *E, real *B,
        const real &n_bg, const real &lambda, const real &guide_field, const real &noise_level, const int *dimX,
        const int *BD, const real *xb, const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_island_coalescence_(real *n_e, real *n_i, real *u_e, real *u_i,
        const real &T0_e, const real &T0_i, real *E, real *B, const real &lambda, const real &eta, const real &psi,
        const real &n_bg, const real &guide_field, const real &noise_level, const int *dimX, const int *BD,
        const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_cssi_(real *n_e, real *n_i, real *u_e, real *u_i, real *T_e, real *T_i, real *E, real *B,
        const real &n_bg, const real &T_bg, const real &delta_e, const real &delta_i,
        const real &B0_e, const real &B0_i, const real &V0_e, const real &V0_i, const real &noise_level,
        const int *dimX, const int *BD, const real *xb, const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_cssi_harris_(real *n_e, real *n_i, real *u_e, real *u_i, real *n_e_bg, real *n_i_bg,
        real *u_e_bg, real *u_i_bg, real &T0_e, real &T0_i, real *E, real *B,
        const real &n_harris, const real &n_bg, const real &lambda, const real &noise_level, const int *dimX, const int *BD,
        const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_lhdi_(real *n_e, real *n_i, real *u_e, real *u_i, real &T0_e, real &T0_i, real *E, real *B,
        const real &n_bg, const real &lambda, const real &noise_level, const int *dimX, const int *BD,
        const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_khi_(real *n_e, real *n_i, real *u_e, real *u_i, real *T_e, real *T_i,
        real *E, real *B, const real &lambda, const real &T0e_T0i, const real &wpe_wce, const real &m_i,
        const int *dimX, const int *BD, const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_whistler_wave_(real *n_e, real *n_i, real *u_e, real *u_i, real *E, real *B,
        const int *dimX, const int *BD, const real *dx, const real *xb_loc);

extern "C" void setup_electromagnetic_vacuum_wave_(real *n_e, real *n_i, real *u_e, real *u_i, real *E, real *B,
        const int *dimX, const int *BD, const real *xb, const real *xe, const real *dx, const real *xb_loc);

extern "C" void setup_electromagnetic_vacuum_wave_3d_plane_polarized_(real *n_e, real *n_i, real *u_e, real *u_i, real *E, real *B,
        const int *dimX, const int *BD, const real *xb, const real *xe, const real *dx, const real *xb_loc);

// wrapper for fortran functions
void init(model::Model *plasma_model[], model::Model *field_solver, const Parameter &p);

} // namespace setup

#endif // SETUP_H_

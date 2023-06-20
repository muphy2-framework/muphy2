/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "framework/definitions.h"

struct Parameter;
struct Species;
struct Block_id;
namespace model {
class Model;
}

namespace output {

// fortran functions
extern "C" void output_extract_distribution_function_slice_(real* f_slice, const real* f, const int* dimX,
                    const int* dimV, const int& max_dimV, const int* BD, const int& direction, const int* slice_coordinates);
extern "C" void output_calc_statistics_fluid_(const real* n, const real* u, const real* T, const real& m,
                    real& mass, real& kinetic_energy, real& thermal_energy, real& momentum_x, real& momentum_y, real& momentum_z,
                    const int* dimX, const int* BD, const real* dx, const int& dimensionality_x, const int& dimensionality_v);
extern "C" void output_calc_statistics_electromagnetic_(const real* E, const real* B,
                    const real& eps0, const real& mu0, real& electric_energy, real& magnetic_energy,
                    real& em_momentum_x, real& em_momentum_y, real& em_momentum_z,
                    const int* dimX, const int* BD, const real* dx, const int& dimensionality_x);

// wrappers for fortran functions
void extract_distribution_function_slice(real* f_slice, const real* f, int max_res_v, const Species &s, const Parameter &p);
void calc_statistics_fluid(const real* n, const real* u, const real* T,
        real& mass, real& kinetic_energy, real& thermal_energy,
        real& momentum_x, real& momentum_y, real& momentum_z, const Species &s, const Parameter &p);
void calc_statistics_electromagnetic(const real* E, const real* B, real& electric_energy, real& magnetic_energy, 
                    real& em_momentum_x, real& em_momentum_y, real& em_momentum_z, const Parameter &p);

// cpp functions
void prepare_and_write(model::Model* plasma_model[], model::Model* field_solver,
        int nplasma_models, real physical_time, int output_number, bool output_vtk,
        const Block_id& block_id, const Parameter &p, int electron_substep_counter=0);

void write_fluid_quantities(real* n[], real* u[], real* T[], real* P[], real* Q_h[],
        real* j, real* E, real* B, real* scheme_id, real* compute_unit, real physical_time,
        int output_number, const Block_id& block_id, const Parameter &p);

void write_distribution_function(model::Model* plasma_model[], real physical_time,
        int output_number, const Block_id& block_id, const Parameter &p);

void write_velocity_distribution(model::Model* plasma_model[], int output_number,
                                 const Block_id& block_id, const Parameter &p);

void set_phase_space_output_parameters(real* xb_loc, real* xb, real* dx, int* coords,
        int* res, int* res_total, int* bd, const Block_id& block_id, const Parameter &p);

bool this_block_phase_space_slice_output(const Block_id& block_id, const Parameter &p);

void write_csv(real* n[], real* u[], real* T[], real* E, real* B, real physical_time,
        int output_number, const Block_id& block_id, const Parameter &p);

void write_restart(model::Model* plasma_model[], model::Model* electromagnetic_model,
                   int nplasma_models, real physical_time, const int output_number,
                   const Block_id& block_id, const Parameter &p, int electron_substep_counter);
}

#endif // OUTPUT_H_

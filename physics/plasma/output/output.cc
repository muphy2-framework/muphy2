/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "output.h"

#include <iostream>

#include "framework/output_write.h"
#include "framework/parameter.h"
#include "framework/block.h"
#include "framework/model.h"
#include "framework/interpolate.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/models/fluid5.h"
#include "physics/plasma/models/mhd.h"
#include "physics/plasma/models/maxwell.h"
#include "physics/plasma/models/poisson.h"
#include "physics/plasma/models/ohm.h"
#include "physics/plasma/converters/converters.h"

namespace output {

void prepare_and_write(model::Model *plasma_model[], model::Model *electromagnetic_model,
                       int nplasma_models, real physical_time, int output_number, bool output_vtk,
                       const Block_id& block_id, const Parameter &p, int electron_substep_counter) {
  real *n[p.nspecies];
  real *u[p.nspecies];
  real *T[p.nspecies];
  real *P[p.nspecies];
  real *Q_h[p.nspecies];

  real *j = new real[3 * p.ncells_x]();
  real *E = new real[3 * p.ncells_x]();
  real *B = new real[3 * p.ncells_x]();
  real *scheme_id = new real[p.ncells_x]();
  real *compute_unit = new real[p.ncells_x]();

  // prepare
  for (int i = 0; i < p.nspecies; ++i) {
    n[i] = new real[p.ncells_x]();
    u[i] = new real[3 * p.ncells_x]();
    T[i] = new real[p.ncells_x]();
    P[i] = new real[6 * p.ncells_x]();
    if (p.heatflux_output) Q_h[i] = new real[10 * p.ncells_x]();

    switch (plasma_model[i]->get_model_id()) {
    case Parameter::kVlasov: {
      model::Vlasov * vlasov = (model::Vlasov *) plasma_model[i];
      #pragma acc update host(vlasov->f[0:p.species[i].ncells_xv])
      convert::f_to_output(vlasov->f, n[i], u[i], T[i], P[i], p.species[i], p);
      if (p.heatflux_output) convert::f_to_heatflux(Q_h[i], vlasov->f, u[i], p.species[i], p);
      if (i==p.nspecies-1) convert::u_to_j(j, u[Parameter::kElectron], u[Parameter::kIon], n[Parameter::kElectron], n[Parameter::kIon], p);
    } break;
    case Parameter::kFluid5: {
      model::Fluid5 * fluid5 = (model::Fluid5 *) plasma_model[i];
      #pragma acc update host(fluid5->five_moments[0:5*p.ncells_x])
      convert::five_moments_to_output(fluid5->five_moments, n[i],
                                      u[i], T[i], P[i], p.species[i].m, p);
      if (i==p.nspecies-1) convert::u_to_j(j, u[Parameter::kElectron], u[Parameter::kIon], n[Parameter::kElectron], n[Parameter::kIon], p);
    } break;
    case Parameter::kFluid10: {
      model::Fluid10 * fluid10 = (model::Fluid10 *) plasma_model[i];
      #pragma acc update host(fluid10->ten_moments[0:10*p.ncells_x])
      convert::ten_moments_to_output(fluid10->ten_moments, n[i], u[i],
                                     T[i], P[i], p.species[i].m, p);
      if (i==p.nspecies-1) convert::u_to_j(j, u[Parameter::kElectron], u[Parameter::kIon], n[Parameter::kElectron], n[Parameter::kIon], p);
    } break;
    case Parameter::kMHD1Temperature: {
      model::MHD1Temperature * mhd = (model::MHD1Temperature *) plasma_model[i];  // TODO: Refactor Output
      convert::five_moments_to_output(mhd->five_moments, n[i],  // TODO, separately for electrons, ions
                                      u[i], T[i], P[i], p.species[Parameter::kIon].m, p);
    } break;
      // other models
    }
  }

  switch (electromagnetic_model->get_model_id()) {
  case Parameter::kMaxwell: {
    model::Maxwell * maxwell = (model::Maxwell *)electromagnetic_model;
    #pragma acc update host(maxwell->E[0:3*p.ncells_x], \
                            maxwell->B[0:3*p.ncells_x])
    for (int i = 0; i < 3 * p.ncells_x; ++i) {
      E[i] = maxwell->E[i];
      B[i] = maxwell->B[i];
    }
  } break;
  case Parameter::kPoisson: {
    model::Poisson * poisson = (model::Poisson *)electromagnetic_model;
    #pragma acc update host(poisson->E[0:3*p.ncells_x], \
                            poisson->B[0:3*p.ncells_x])
    for (int i = 0; i < 3 * p.ncells_x; ++i) {
      E[i] = poisson->E[i];
      B[i] = poisson->B[i];
    }
  } break;
  case Parameter::kOhm: {
    model::Ohm * ohm = (model::Ohm *)electromagnetic_model;
    for (int i = 0; i < 3 * p.ncells_x; ++i) {
      E[i] = ohm->E[i];
      B[i] = ohm->B[i];
      j[i] = ohm->J[i];
    }
  } break;
    // other field solvers
  }

  // output
  write_csv(n, u, T, E, B, physical_time, output_number, block_id, p);

  if (output_vtk) {

    write_fluid_quantities(n, u, T, P, Q_h, j, E, B, scheme_id, compute_unit,
                           physical_time, output_number, block_id, p);

    if (p.phase_space_output == "slice" && this_block_phase_space_slice_output(block_id, p)) {
      write_distribution_function(plasma_model, physical_time,
                                  output_number, block_id, p);
    }

    if (p.restart_output || p.phase_space_output == "full" || p.phase_space_output == "quarter") {
      // restart output is repurposed also for full distribution function output
      write_restart(plasma_model, electromagnetic_model, nplasma_models, physical_time, output_number,
                    block_id, p, electron_substep_counter);
    }
  }

  // clean up
  delete[] B;
  delete[] E;
  delete[] j;
  delete[] compute_unit;
  delete[] scheme_id;
  for (int i = p.nspecies-1; i >= 0; --i) {
    delete[] P[i];
    delete[] T[i];
    delete[] u[i];
    delete[] n[i];
  }
}


void write_fluid_quantities(real* n[], real* u[], real* T[], real* P[], real* Q_h[],
        real* j, real* E, real* B, real* scheme_id, real* compute_unit, real physical_time,
        int output_number, const Block_id& block_id, const Parameter &p) {

  std::string field_names_species_dependent[] =
      {"n", "ux", "uy", "uz", "T", "Pxx", "Pxy", "Pxz", "Pyy", "Pyz", "Pzz"};
  std::string field_names_heatflux[] =
      {"Qxxx", "Qxxy", "Qxxz", "Qxyy", "Qxyz", "Qxzz", "Qyyy", "Qyyz", "Qyzz", "Qzzz"};
  std::string field_names_other[] =
      {"jx", "jy", "jz", "Ex", "Ey", "Ez", "Bx", "By", "Bz", "scheme_id", "compute_unit"};

  int nfields_species_dependent = sizeof(field_names_species_dependent) /
                                  sizeof(field_names_species_dependent[0]);
  int nfields_heatflux = 0;
  if (p.heatflux_output) {
    nfields_heatflux = sizeof(field_names_heatflux) / sizeof(field_names_heatflux[0]);
  }
  int nfields_other = sizeof(field_names_other) / sizeof(field_names_other[0]);
  int nfields_total = (nfields_species_dependent + nfields_heatflux) * p.nspecies + nfields_other;
  std::string field_names[nfields_total];
  int i = 0;
  for (int ifield = 0; ifield < nfields_species_dependent; ++ifield) {
    for (int ispecies = 0; ispecies < p.nspecies; ++ispecies) {
      field_names[i] = field_names_species_dependent[ifield] + "_" +
                       p.species_names_str[ispecies];
      ++i;
    }
  }
  for (int ifield = 0; ifield < nfields_heatflux; ++ifield) {
    for (int ispecies = 0; ispecies < p.nspecies; ++ispecies) {
      field_names[i] = field_names_heatflux[ifield] + "_" +
                       p.species_names_str[ispecies];
      ++i;
    }
  }
  for (int ifield = 0; ifield < nfields_other; ++ifield) {
    field_names[(nfields_species_dependent + nfields_heatflux) * p.nspecies + ifield] =
        field_names_other[ifield];
  }

  const real *fields[nfields_total];

  int ifields = 0;
  for (int i = 0; i < p.nspecies; ++i) {
    fields[ifields] = &n[i][0];
    ++ifields;
  }
  for (int k = 0; k < 3; ++k) {
    for (int i = 0; i < p.nspecies; ++i) {
      fields[ifields] = &u[i][k * p.ncells_x];
      ++ifields;
    }
  }
  for (int i = 0; i < p.nspecies; ++i) {
    fields[ifields] = &T[i][0];
    ++ifields;
  }
  for (int k = 0; k < 6; ++k) {
    for (int i = 0; i < p.nspecies; ++i) {
      fields[ifields] = &P[i][k * p.ncells_x];
      ++ifields;
    }
  }
  if (p.heatflux_output) {
    for (int k = 0; k < 10; ++k) {
      for (int i = 0; i < p.nspecies; ++i) {
        fields[ifields] = &Q_h[i][k * p.ncells_x];
        ++ifields;
      }
    }
  }
  for (int k = 0; k < 3; ++k) {
    fields[ifields] = &j[k * p.ncells_x];
    ++ifields;
  }
  for (int k = 0; k < 3; ++k) {
    fields[ifields] = &E[k * p.ncells_x];
    ++ifields;
  }
  for (int k = 0; k < 3; ++k) {
    fields[ifields] = &B[k * p.ncells_x];
    ++ifields;
  }
  for (int cell = 0; cell < p.ncells_x; ++cell) {
    scheme_id[cell] = block_id.my_scheme_id;
  }
  fields[ifields] = scheme_id;
  ++ifields;
  if (p.nslots_per_node > 0) {
    for (int cell = 0; cell < p.ncells_x; ++cell) {
      compute_unit[cell] = block_id.compute_node*p.nslots_per_node + block_id.selected_device;
    }
  }
  fields[ifields] = compute_unit;

  write_vti_file(fields, nfields_total, field_names, output_number,
                 p.output_directory, p.xb_loc, p.dx, block_id.my_coords, p.res_x, p.bd);

  if (block_id.my_coords[0] == 0 && block_id.my_coords[1] == 0 && block_id.my_coords[2] == 0) {
    write_pvti_file(nfields_total, field_names, physical_time, output_number,
                    p.output_directory, p.xb, p.dx, p.res_x, p.res_x_total);
  }
}


void write_distribution_function(model::Model* plasma_model[], real physical_time,
        int output_number, const Block_id& block_id, const Parameter &p) {

  int coords[3], res[3], res_total[3], bd[3];
  real xb_loc[3], xb[3], dx[3];
  set_phase_space_output_parameters(xb_loc, xb, dx, coords, res, res_total, bd, block_id, p);

  std::string field_names[p.nspecies];
  const real *fields_f[p.nspecies];
  real *f_slice[p.nspecies];

  for (int i = 0; i < p.nspecies; ++i) {
    f_slice[i] = new real[(res[0]+2*bd[0])*res[1]]();
    if (plasma_model[i]->get_model_id() == Parameter::kVlasov) {
      model::Vlasov * vlasov = (model::Vlasov *) plasma_model[i];
      #pragma acc update host(vlasov->f[0:p.species[i].ncells_xv])
      extract_distribution_function_slice(f_slice[i], vlasov->f, res[1], p.species[i], p);
    }
    fields_f[i] = &f_slice[i][0];

    field_names[i] = "f_" + p.species_names_str[i];
  }

  write_vti_file(fields_f, p.nspecies, field_names, output_number,
                 p.output_directory + "/phase_space/", xb_loc, dx, coords, res, bd);

  if (block_id.my_coords[0] == 0 && block_id.my_coords[1] == 0 && block_id.my_coords[2] == 0) {
    write_pvti_file(p.nspecies, field_names, physical_time, output_number,
                    p.output_directory + "/phase_space/", xb, dx, res, res_total);
  }

  for (int i = p.nspecies-1; i <= 0; --i) {
    delete[] f_slice[i];
  }
}


void write_csv(real* n[], real* u[], real* T[], real* E, real* B, real physical_time,
        int output_number, const Block_id& block_id, const Parameter &p) {

  std::string field_names_species_dependent[] =
      {"mass", "kinetic energy", "thermal energy", "momentum_x", "momentum_y", "momentum_z"};
  std::string field_names_other[] =
      {"time", "electric energy", "magnetic energy", "em_momentum_x", "em_momentum_y", "em_momentum_z"};

  int nfields_species_dependent = sizeof(field_names_species_dependent) /
                                  sizeof(field_names_species_dependent[0]);
  int nfields_other = sizeof(field_names_other) / sizeof(field_names_other[0]);
  int nfields_total = nfields_species_dependent * p.nspecies + nfields_other;

  std::string field_names[nfields_total];
  for (int ifield = 0; ifield < nfields_other; ++ifield) {
    field_names[ifield] = field_names_other[ifield];
  }
  for (int ispecies = 0; ispecies < p.nspecies; ++ispecies) {
    for (int ifield = 0; ifield < nfields_species_dependent; ++ifield) {
      field_names[nfields_other+nfields_species_dependent*ispecies+ifield] = 
          field_names_species_dependent[ifield] + " " + p.species_names_str[ispecies];
    }
  }

  real fields_text_output[nfields_total];

  fields_text_output[0] = physical_time;

  real *E_cc = new real[3 * p.ncells_x]();
  real *B_cc = new real[3 * p.ncells_x]();
  ipol::yee_face_to_centered(E, E_cc, p);
  ipol::yee_edge_to_centered(B, B_cc, p);

  calc_statistics_electromagnetic(E_cc, B_cc, fields_text_output[1], fields_text_output[2],
              fields_text_output[3], fields_text_output[4], fields_text_output[5], p);

  delete[] B_cc;
  delete[] E_cc;

  for (int ispecies = 0; ispecies < p.nspecies; ++ispecies) {
    calc_statistics_fluid(n[ispecies], u[ispecies], T[ispecies],
            fields_text_output[6+ispecies*nfields_species_dependent], fields_text_output[7+ispecies*nfields_species_dependent],
            fields_text_output[8+ispecies*nfields_species_dependent], fields_text_output[9+ispecies*nfields_species_dependent],
            fields_text_output[10+ispecies*nfields_species_dependent], fields_text_output[11+ispecies*nfields_species_dependent],
            p.species[ispecies], p);
  }

  write_text_file(fields_text_output, nfields_total, field_names, output_number,
                   p.output_directory, block_id.my_coords);
}


void write_restart(model::Model* plasma_model[], model::Model* electromagnetic_model,
                   int nplasma_models, real physical_time, const int output_number,
                   const Block_id& block_id, const Parameter &p, int electron_substep_counter) {

  int nfields = nplasma_models,
      nscalars = 0;

  if (electromagnetic_model->get_model_id() == Parameter::kMaxwell) {
    nscalars += 2; // dt_maxwell, substeps
    nfields  += 3; // E, B, B_ipol
  } else if (electromagnetic_model->get_model_id() == Parameter::kPoisson) {
    nfields  += 1; // phi
  }

  if (electron_substep_counter != 0) {
    nscalars++;
  }

  // fill scalars array
  real scalars[nscalars];
  if (electromagnetic_model->get_model_id() == Parameter::kMaxwell) {
    scalars[0] = ((model::Maxwell *)electromagnetic_model)->dt;
    scalars[1] = ((model::Maxwell *)electromagnetic_model)->substeps;
  }

  if (electron_substep_counter != 0) {
    scalars[2] = electron_substep_counter;
  }

  // fill fields array
  const real* fields[nfields];
  int field_sizes[nfields];

  for (int i = 0; i < nplasma_models; ++i) {
    switch (plasma_model[i]->get_model_id()) {
    case Parameter::kVlasov: {
      if (i >= p.nspecies) {
        std::cerr<<std::endl<<"Error in output.cc: Please check order of passed models."<<std::endl;
        std::exit(EXIT_FAILURE);
      }
      fields[i] = &(((model::Vlasov *)plasma_model[i])->f[0]);
      field_sizes[i] = p.species[i].ncells_xv;
    } break;
    case Parameter::kFluid10: {
#ifdef _OPENACC
      if (i >= p.nspecies) { // otherwise already updated
        model::Fluid10 * fluid10 = (model::Fluid10 *) plasma_model[i];
        #pragma acc update host(fluid10->ten_moments[0:10*p.ncells_x])
      }
#endif
      fields[i] = &(((model::Fluid10 *)plasma_model[i])->ten_moments[0]);
      field_sizes[i] = 10*p.ncells_x;
    } break;
    case Parameter::kFluid5: {
      fields[i] = &(((model::Fluid5 *)plasma_model[i])->five_moments[0]);
      field_sizes[i] = 5*p.ncells_x;
    } break;
      // other models
    }
  }

  if (electromagnetic_model->get_model_id() == Parameter::kMaxwell) {
    fields[nplasma_models] = &(((model::Maxwell *)electromagnetic_model)->E[0]);
    field_sizes[nplasma_models] = 3*p.ncells_x;
    fields[nplasma_models+1] = &(((model::Maxwell *)electromagnetic_model)->B[0]);
    field_sizes[nplasma_models+1] = 3*p.ncells_x;
    fields[nplasma_models+2] = &(((model::Maxwell *)electromagnetic_model)->B_ipol[0]);
    field_sizes[nplasma_models+2] = 3*p.ncells_x;
  } else if (electromagnetic_model->get_model_id() == Parameter::kPoisson) {
    fields[nplasma_models] = &(((model::Poisson *)electromagnetic_model)->phi[0]);
    field_sizes[nplasma_models] = p.ncells_x;
  }

  // output
  if (p.phase_space_output == "full" || (p.phase_space_output == "quarter" && 
        block_id.my_coords[0] < p.nproc[0]-p.nproc[0]/2 &&
        block_id.my_coords[1] >= p.nproc[1]/2 &&
        block_id.my_coords[2] == p.nproc[2]/2)) {
    write_restart_file(block_id.my_scheme_id, physical_time,
                       scalars, nscalars, fields, field_sizes, nfields,
                       p.output_directory, block_id.my_coords, output_number);
  } else if (p.restart_output) {
    write_restart_file(block_id.my_scheme_id, physical_time,
                       scalars, nscalars, fields, field_sizes, nfields,
                       p.output_directory, block_id.my_coords);
  }
}


bool this_block_phase_space_slice_output(const Block_id& block_id, const Parameter &p) {

  int dir = p.phase_space_output_direction;
  int index_x = block_id.my_coords[0]*p.res_x[0],
      index_y = block_id.my_coords[1]*p.res_x[1],
      index_z = block_id.my_coords[2]*p.res_x[2];

  // check if slice lies in the respective block, otherwise return false
  if (dir == 0) { // x
    if (p.phase_space_output_slice_index[0] < index_y ||
        p.phase_space_output_slice_index[0] >= index_y+p.res_x[1] ||
        p.phase_space_output_slice_index[1] < index_z ||
        p.phase_space_output_slice_index[1] >= index_z+p.res_x[2]) {
        return false;
    }
  } else if (dir == 1) { // y
    if (p.phase_space_output_slice_index[0] < index_x ||
        p.phase_space_output_slice_index[0] >= index_x+p.res_x[0] ||
        p.phase_space_output_slice_index[1] < index_z ||
        p.phase_space_output_slice_index[1] >= index_z+p.res_x[2]) {
        return false;
    }
  } else if (dir == 2) { // z
    if (p.phase_space_output_slice_index[0] < index_x ||
        p.phase_space_output_slice_index[0] >= index_x+p.res_x[0] ||
        p.phase_space_output_slice_index[1] < index_y ||
        p.phase_space_output_slice_index[1] >= index_y+p.res_x[1]) {
        return false;
    }
  } else { // no legit direction
    return false;
  }

  return true;
}


void set_phase_space_output_parameters(real* xb_loc, real* xb, real* dx, int* coords,
                    int* res, int* res_total, int* bd, const Block_id& block_id, const Parameter &p) {

  int dir = p.phase_space_output_direction;

  xb_loc[0] = p.res_x[dir]*block_id.my_coords[dir]; xb_loc[1] = 0; xb_loc[2] = 0;
  xb[0] = 0; xb[1] = 0; xb[2] = 0;
  dx[0] = 1; dx[1] = 1; dx[2] = 1;
  res[0] = p.res_x[dir];
  res[1] = 0;
  for (int i = 0; i < p.nspecies; ++i) { // find maximum velocity resolution
    if (p.species[i].res_v[dir] > res[1]) {
      res[1] = p.species[i].res_v[dir];
    }
  }
  res[2] = 1;
  res_total[0] = p.res_x_total[dir]; res_total[1] = res[1]; res_total[2] = 1;
  bd[0] = p.bd[dir]; bd[1] = 0; bd[2] = 0;
  coords[0] = block_id.my_coords[dir]; coords[1] = 0; coords[2] = 0;
}


// fortran wrappers
void extract_distribution_function_slice(real* f_slice, const real* f, int max_res_v, const Species &s, const Parameter &p) {

  output_extract_distribution_function_slice_(f_slice, f, p.res_x_minus_one, s.res_v, max_res_v, p.bd,
                    p.phase_space_output_direction+1, p.phase_space_output_slice_index);
}


void calc_statistics_fluid(const real* n, const real* u, const real* T,
        real& mass, real& kinetic_energy, real& thermal_energy,
        real& momentum_x, real& momentum_y, real& momentum_z, const Species &s, const Parameter &p) {

  output_calc_statistics_fluid_(n, u, T, s.m, mass, kinetic_energy, thermal_energy,
                                momentum_x, momentum_y, momentum_z, p.res_x_minus_one,
                                p.bd, p.dx, p.dimensionality_x, p.dimensionality_v);
}


void calc_statistics_electromagnetic(const real* E, const real* B, real& electric_energy, real& magnetic_energy, 
                    real& em_momentum_x, real& em_momentum_y, real& em_momentum_z, const Parameter &p) {

  output_calc_statistics_electromagnetic_(E, B, p.eps0, p.mu0, electric_energy, magnetic_energy,
                                          em_momentum_x, em_momentum_y, em_momentum_z, p.res_x_minus_one,
                                          p.bd, p.dx, p.dimensionality_x);
}

} // namespace output

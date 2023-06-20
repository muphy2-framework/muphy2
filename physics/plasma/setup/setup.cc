/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "setup.h"

#include <iostream>

#include "framework/utilities.h"
#include "framework/parameter.h"
#include "framework/model.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/models/fluid5.h"
#include "physics/plasma/models/mhd.h"
#include "physics/plasma/models/maxwell.h"
#include "physics/plasma/models/poisson.h"
#include "physics/plasma/models/ohm.h"
#include "physics/plasma/converters/converters.h"

namespace setup {

void init(model::Model *plasma_model[], model::Model *electromagnetic_model, const Parameter &p) {

  // temporary fields for initialization
  real* f[p.nspecies]; // not allocated if fluid initialization, TODO: do not allocate for Vlasov model init
  real* n[p.nspecies];
  real* u[p.nspecies];
  real* T[p.nspecies];

  real* n_bg[p.nspecies]; // separate background values, only allocated when needed
  real* u_bg[p.nspecies];
  real* T_bg[p.nspecies];

  for (int i = 0; i < p.nspecies; ++i) {
    n[i] = new real[p.ncells_x]();
    u[i] = new real[3*p.ncells_x]();
    T[i] = new real[p.ncells_x];
    util::array_fill(T[i], p.species[i].T0, p.ncells_x); // default temperature init
  }

  real* E = new real[3*p.ncells_x]();
  real* B = new real[3*p.ncells_x]();

  // distribution function setups two species
  bool setup_type_distribution_function_two_species = true;

  if (p.setup == "landau_damping") {
    const std::string necessary_keys[1] = {"alpha"};
    util::check_umap_keys(necessary_keys, 1, p);
    for (int i = 0; i < p.nspecies; ++i) {
      f[i] = new real[p.species[i].ncells_xv]();
    }
    setup_landau_damping_(f[Parameter::kElectron], f[Parameter::kIon], E, B,
                        p.species[Parameter::kElectron].m, p.species[Parameter::kElectron].T0,
                        p.eps0, p.setup_var.at("alpha"), p.res_x_minus_one,
                        p.species[Parameter::kElectron].res_v, p.species[Parameter::kIon].res_v, p.bd, p.dx,
                        p.species[Parameter::kElectron].vb, p.species[Parameter::kElectron].dv,
                        p.species[Parameter::kIon].dv, p.dimensionality_x, p.dimensionality_v, p.xb_loc);
  } else if (p.setup == "two_stream_instability") {
    const std::string necessary_keys[1] = {"alpha"};
    util::check_umap_keys(necessary_keys, 1, p);
    for (int i = 0; i < p.nspecies; ++i) {
      f[i] = new real[p.species[i].ncells_xv]();
    }
    setup_two_stream_instability_(f[Parameter::kElectron], f[Parameter::kIon], E, B,
                        p.species[Parameter::kElectron].m, p.species[Parameter::kElectron].T0,
                        p.setup_var.at("alpha"), p.res_x_minus_one,
                        p.species[Parameter::kElectron].res_v, p.species[Parameter::kIon].res_v, p.bd, p.dx,
                        p.species[Parameter::kElectron].vb, p.species[Parameter::kElectron].dv,
                        p.species[Parameter::kIon].dv, p.dimensionality_v, p.xb_loc);
  } else {
    setup_type_distribution_function_two_species = false;
  }
    

  // five moment setups two species
  bool setup_type_five_moments_two_species = true;

  if (p.setup == "orszag_tang") {
    const std::string necessary_keys[3] = {"L","n_bg","delta_u"};
    util::check_umap_keys(necessary_keys, 3, p);
    setup_orszag_tang_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        E, B, p.c0, p.setup_var.at("L"), p.setup_var.at("n_bg"), p.setup_var.at("delta_u"),
                        p.res_x_minus_one, p.bd, p.dx, p.xb_loc);
  } else if (p.setup == "orszag_tang_3d") {
    const std::string necessary_keys[5] = {"L","L_z","n_bg","delta_u","cross_helicity"};
    util::check_umap_keys(necessary_keys, 5, p);
    setup_orszag_tang_3d_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        E, B, p.setup_var.at("L"), p.setup_var.at("L_z"), p.setup_var.at("n_bg"), p.setup_var.at("delta_u"),
                        p.setup_var.at("cross_helicity"), p.res_x_minus_one, p.bd, p.dx, p.xb_loc);
  } else if (p.setup == "double_harris_sheet") {
    const std::string necessary_keys[5] = {"lambda","psi","n_bg","guide_field","noise_level"};
    util::check_umap_keys(necessary_keys, 5, p);
    setup_double_harris_sheet_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        p.species[Parameter::kElectron].T0, p.species[Parameter::kIon].T0, E, B,
                        p.setup_var.at("lambda"), p.setup_var.at("psi"), p.setup_var.at("n_bg"), p.setup_var.at("guide_field"),
                        p.setup_var.at("noise_level"), p.res_x_minus_one, p.bd, p.xb, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "island_coalescence") {
    const std::string necessary_keys[6] = {"lambda","eta","psi","n_bg","guide_field","noise_level"};
    util::check_umap_keys(necessary_keys, 6, p);
    setup_island_coalescence_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        p.species[Parameter::kElectron].T0, p.species[Parameter::kIon].T0, E, B,
                        p.setup_var.at("lambda"), p.setup_var.at("eta"), p.setup_var.at("psi"),
                        p.setup_var.at("n_bg"), p.setup_var.at("guide_field"), p.setup_var.at("noise_level"),
                        p.res_x_minus_one, p.bd, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "cssi") {
    const std::string necessary_keys[9] = {"n_bg","T_bg","delta_e","delta_i","B0_e","B0_i","V0_e","V0_i","noise_level"};
    util::check_umap_keys(necessary_keys, 9, p);
    setup_cssi_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        T[Parameter::kElectron], T[Parameter::kIon], E, B, p.setup_var.at("n_bg"), p.setup_var.at("T_bg"),
                        p.setup_var.at("delta_e"), p.setup_var.at("delta_i"), p.setup_var.at("B0_e"), p.setup_var.at("B0_i"),
                        p.setup_var.at("V0_e"), p.setup_var.at("V0_i"), p.setup_var.at("noise_level"),
                        p.res_x_minus_one, p.bd, p.xb, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "lhdi") {
    const std::string necessary_keys[3] = {"n_bg","lambda","noise_level"};
    util::check_umap_keys(necessary_keys, 3, p);
    setup_lhdi_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        p.species[Parameter::kElectron].T0, p.species[Parameter::kIon].T0, E, B,
                        p.setup_var.at("n_bg"), p.setup_var.at("lambda"), p.setup_var.at("noise_level"),
                        p.res_x_minus_one, p.bd, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "khi") {
    const std::string necessary_keys[3] = {"lambda","T0e_T0i","wpe_wce"};
    util::check_umap_keys(necessary_keys, 3, p);
    setup_khi_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        T[Parameter::kElectron], T[Parameter::kIon], E, B, p.setup_var.at("lambda"),
                        p.setup_var.at("T0e_T0i"), p.setup_var.at("wpe_wce"), p.species[Parameter::kIon].m,
                        p.res_x_minus_one, p.bd, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "whistler_wave") {
    setup_whistler_wave_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        E, B, p.res_x_minus_one, p.bd, p.dx, p.xb_loc);
  } else if (p.setup == "electromagnetic_vacuum_wave") {
    setup_electromagnetic_vacuum_wave_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        E, B, p.res_x_minus_one, p.bd, p.xb, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "electromagnetic_vacuum_wave_3d_plane_polarized") {
    setup_electromagnetic_vacuum_wave_3d_plane_polarized_(n[Parameter::kElectron], n[Parameter::kIon],
                        u[Parameter::kElectron], u[Parameter::kIon], E, B, p.res_x_minus_one, p.bd, p.xb, p.xe, p.dx, p.xb_loc);
  } else {
    setup_type_five_moments_two_species = false;
  }


  // five moment setups two species with separate background plasma
  bool setup_type_five_moments_two_species_bg = true;

  if (p.setup == "harris_sheet") {
    const std::string necessary_keys[9] = {"T_bg_e","T_bg_i","lambda","psi","n_bg","guide_field",
                                           "noise_level","drifting_background","sine_perturbation"};
    util::check_umap_keys(necessary_keys, 9, p);
    for (int i = 0; i < p.nspecies; ++i) {
      n_bg[i] = new real[p.ncells_x]();
      u_bg[i] = new real[3*p.ncells_x]();
      T_bg[i] = new real[p.ncells_x];
    }
    util::array_fill(T_bg[Parameter::kElectron], p.setup_var.at("T_bg_e"), p.ncells_x);
    util::array_fill(T_bg[Parameter::kIon], p.setup_var.at("T_bg_i"), p.ncells_x);
    setup_harris_sheet_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        n_bg[Parameter::kElectron], n_bg[Parameter::kIon], u_bg[Parameter::kElectron], u_bg[Parameter::kIon],
                        p.species[Parameter::kElectron].T0, p.species[Parameter::kIon].T0, E, B,
                        p.setup_var.at("lambda"), p.setup_var.at("psi"), p.setup_var.at("n_bg"),
                        p.setup_var.at("guide_field"), p.setup_var.at("noise_level"),p.setup_var_bool.at("drifting_background"),
                        p.setup_var_bool.at("sine_perturbation"), p.res_x_minus_one, p.bd, p.xb, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "double_harris_sheet_gaussian_hump") {
    const std::string necessary_keys[4] = {"n_bg","lambda","guide_field","noise_level"};
    util::check_umap_keys(necessary_keys, 4, p);
    for (int i = 0; i < p.nspecies; ++i) {
      n_bg[i] = new real[p.ncells_x]();
      u_bg[i] = new real[3*p.ncells_x]();
      T_bg[i] = new real[p.ncells_x];
    }
    util::array_fill(T_bg[Parameter::kElectron], p.species[Parameter::kElectron].T0, p.ncells_x);
    util::array_fill(T_bg[Parameter::kIon], p.species[Parameter::kIon].T0, p.ncells_x);
    setup_double_harris_sheet_gaussian_hump_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                        n_bg[Parameter::kElectron], n_bg[Parameter::kIon], u_bg[Parameter::kElectron], u_bg[Parameter::kIon],
                        p.species[Parameter::kElectron].T0, p.species[Parameter::kIon].T0, E, B,
                        p.setup_var.at("n_bg"), p.setup_var.at("lambda"), p.setup_var.at("guide_field"),
                        p.setup_var.at("noise_level"), p.res_x_minus_one, p.bd, p.xb, p.xe, p.dx, p.xb_loc);
  } else if (p.setup == "cssi_harris") {
    const std::string necessary_keys[6] = {"T_bg_e","T_bg_i","lambda","n_harris","n_bg","noise_level"};
    util::check_umap_keys(necessary_keys, 6, p);
    for (int i = 0; i < p.nspecies; ++i) {
      n_bg[i] = new real[p.ncells_x]();
      u_bg[i] = new real[3*p.ncells_x]();
      T_bg[i] = new real[p.ncells_x];
    }
    util::array_fill(T_bg[Parameter::kElectron], p.setup_var.at("T_bg_e"), p.ncells_x);
    util::array_fill(T_bg[Parameter::kIon], p.setup_var.at("T_bg_i"), p.ncells_x);
    setup_cssi_harris_(n[Parameter::kElectron], n[Parameter::kIon], u[Parameter::kElectron], u[Parameter::kIon],
                       n_bg[Parameter::kElectron], n_bg[Parameter::kIon], u_bg[Parameter::kElectron], u_bg[Parameter::kIon],
                       p.species[Parameter::kElectron].T0, p.species[Parameter::kIon].T0, E, B,
                       p.setup_var.at("n_harris"), p.setup_var.at("n_bg"), p.setup_var.at("lambda"),
                       p.setup_var.at("noise_level"), p.res_x_minus_one, p.bd, p.xe, p.dx, p.xb_loc);
  } else {
    setup_type_five_moments_two_species_bg = false;
  }


  // convert distribution function init to plasma models
  if (setup_type_distribution_function_two_species) {
    
    const int array_size[3] = {p.res_x[0]+2*p.bd[0], p.res_x[1]+2*p.bd[1], p.res_x[2]+2*p.bd[2]};
    for (int i = 0; i < p.nspecies; ++i) {
      switch (plasma_model[i]->get_model_id()) {
      case Parameter::kVlasov: {
        model::Vlasov *vlasov = (model::Vlasov *)plasma_model[i];
        for (int k = 0; k < p.species[i].ncells_xv; ++k) {
          vlasov->f[k] = (f[i])[k];
        }
        #pragma acc update device(vlasov->f[0:p.species[i].ncells_xv])
      } break;
      case Parameter::kFluid5: {
        model::Fluid5 *fluid5 = (model::Fluid5 *)plasma_model[i];
        convert::f_to_five_moments(fluid5->five_moments, f[i], array_size, p.species[i], p);
      } break;
      case Parameter::kFluid10: {
        model::Fluid10 *fluid10 = (model::Fluid10 *)plasma_model[i];
        convert::f_to_ten_moments(fluid10->ten_moments, f[i], array_size, p.species[i], p);
      } break;
      // TODO MHD
      }
    }
  }

  // convert fluid quantities init to plasma models
  else if (setup_type_five_moments_two_species) {

    for (int i = 0; i < p.nspecies; ++i) {
     
      switch (plasma_model[i]->get_model_id()) {
      case Parameter::kVlasov: {
        model::Vlasov *vlasov = (model::Vlasov *)plasma_model[i];
        convert::fluid_quantities_to_f(vlasov->f, n[i], u[i], T[i], p.species[i], p);
      } break;
      case Parameter::kFluid5: {
        model::Fluid5 *fluid5 = (model::Fluid5 *)plasma_model[i];
        convert::fluid_quantities_to_five_moments(fluid5->five_moments, n[i], u[i], T[i], p.species[i].m, p);
      } break;
      case Parameter::kFluid10: {
        model::Fluid10 *fluid10 = (model::Fluid10 *)plasma_model[i];
        convert::fluid_quantities_to_ten_moments(fluid10->ten_moments, n[i], u[i], T[i], p.species[i].m, p);
      } break;
      case Parameter::kMHD1Temperature: {
        model::MHD1Temperature *mhd = (model::MHD1Temperature *)plasma_model[i]; // i=0
        convert::fluid_quantities_to_mhd(mhd->five_moments, n[Parameter::kElectron], n[Parameter::kIon],
                                         u[Parameter::kElectron], u[Parameter::kIon],
                                         T[Parameter::kElectron], T[Parameter::kIon],
                                         p.species[Parameter::kElectron].m, p.species[Parameter::kIon].m, p);
      } break;
      default:
        std::cerr << "Error in setup: Plasma model not found." << std::endl;
      }
    }
  }

  // convert fluid quantities init with separate background to plasma models
  else if (setup_type_five_moments_two_species_bg) {

    for (int i = 0; i < p.nspecies; ++i) {
     
      switch (plasma_model[i]->get_model_id()) {
      case Parameter::kVlasov: {
        model::Vlasov *vlasov = (model::Vlasov *)plasma_model[i];
        convert::fluid_quantities_to_f_separate_background(vlasov->f, n[i], u[i], T[i],
                                            n_bg[i], u_bg[i], T_bg[i], p.species[i], p);
      } break;
      case Parameter::kFluid5: {
        model::Fluid5 *fluid5 = (model::Fluid5 *)plasma_model[i];
        convert::fluid_quantities_to_five_moments(fluid5->five_moments, n[i], u[i], T[i], p.species[i].m, p);
        real* five_moments_bg = new real[5*p.ncells_x]();
        convert::fluid_quantities_to_five_moments(five_moments_bg, n_bg[i], u_bg[i], T_bg[i], p.species[i].m, p);
        util::array_add_equal(fluid5->five_moments, five_moments_bg, 5*p.ncells_x);
        delete[] five_moments_bg;
      } break;
      case Parameter::kFluid10: {
        model::Fluid10 *fluid10 = (model::Fluid10 *)plasma_model[i];
        convert::fluid_quantities_to_ten_moments(fluid10->ten_moments, n[i], u[i], T[i], p.species[i].m, p);
        real* ten_moments_bg = new real[10*p.ncells_x]();
        convert::fluid_quantities_to_ten_moments(ten_moments_bg, n_bg[i], u_bg[i], T_bg[i], p.species[i].m, p);
        util::array_add_equal(fluid10->ten_moments, ten_moments_bg, 10*p.ncells_x);
        delete[] ten_moments_bg;
      } break;
      case Parameter::kMHD1Temperature: {
        model::MHD1Temperature *mhd = (model::MHD1Temperature *)plasma_model[i]; // i=0
        convert::fluid_quantities_to_mhd(mhd->five_moments, n[Parameter::kElectron], n[Parameter::kIon],
                                         u[Parameter::kElectron], u[Parameter::kIon],
                                         T[Parameter::kElectron], T[Parameter::kIon],
                                         p.species[Parameter::kElectron].m, p.species[Parameter::kIon].m, p);
        real* five_moments_bg = new real[5*p.ncells_x]();
        convert::fluid_quantities_to_mhd(mhd->five_moments, n_bg[Parameter::kElectron], n_bg[Parameter::kIon],
                                         u_bg[Parameter::kElectron], u_bg[Parameter::kIon],
                                         T_bg[Parameter::kElectron], T_bg[Parameter::kIon],
                                         p.species[Parameter::kElectron].m, p.species[Parameter::kIon].m, p);
        util::array_add_equal(mhd->five_moments, five_moments_bg, 5*p.ncells_x);
        delete[] five_moments_bg;
      } break;
      default:
        std::cerr << "Error in setup: Plasma model not found." << std::endl;
      }
    }
  }

  // convert electro-magnetic fields init to models
  switch (electromagnetic_model->get_model_id()) {
  case Parameter::kMaxwell: {
    model::Maxwell *maxwell = (model::Maxwell *)electromagnetic_model;
    for (int i = 0; i < 3 * p.ncells_x; ++i) {
      maxwell->E[i] = E[i];
      maxwell->B[i] = B[i];
    }
    #pragma acc update device(maxwell->E[0:3*p.ncells_x], \
                              maxwell->B[0:3*p.ncells_x])
  } break;
  case Parameter::kPoisson: {
    model::Poisson *poisson = (model::Poisson *)electromagnetic_model;
    for (int i = 0; i < 3 * p.ncells_x; ++i) {
      poisson->E[i] = E[i];
      poisson->B[i] = B[i];
    }
    #pragma acc update device(poisson->E[0:3*p.ncells_x], \
                              poisson->B[0:3*p.ncells_x])
  } break;
  case Parameter::kOhm: {
    model::Ohm *ohm = (model::Ohm *)electromagnetic_model;
    for (int i = 0; i < 3 * p.ncells_x; ++i) {
      ohm->E[i] = E[i]; // needed?
      ohm->B[i] = B[i];
    }
  } break;
  default:
    std::cerr << "Error in setup: Electromagnetic model not found." << std::endl;
  }

  delete[] B;
  delete[] E;

  for (int i = 0; i < p.nspecies; ++i) {
    delete[] n[i];
    delete[] u[i];
    delete[] T[i];
  }

  if (setup_type_five_moments_two_species_bg) {
    for (int i = 0; i < p.nspecies; ++i) {
      delete[] n_bg[i];
      delete[] u_bg[i];
      delete[] T_bg[i];
    }
  }

  if (setup_type_distribution_function_two_species) {
     for (int i = 0; i < p.nspecies; ++i) {
      delete[] f[i];
    }
  }

  if (!setup_type_distribution_function_two_species &&
      !setup_type_five_moments_two_species &&
      !setup_type_five_moments_two_species_bg) {
    std::cerr << std::endl << "Error: The specified setup does not exist." << std::endl;
  }

} // void init

} // namespace setup


/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "select_scheme.h"

#include <iostream>

#include "framework/parameter.h"
#include "framework/block.h"
#include "framework/scheme.h"
#include "framework/model.h"
#include "framework/criterion.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/models/fluid5.h"
#include "physics/plasma/models/mhd.h"
#include "physics/plasma/models/maxwell.h"
#include "physics/plasma/models/ohm.h"
#include "physics/plasma/schemes/VeViM.h"
#include "physics/plasma/schemes/VeViMdeltaf.h"
#include "physics/plasma/schemes/VeViP.h"
#include "physics/plasma/schemes/VeViPdeltaf.h"
#include "physics/plasma/schemes/F10eViM.h"
#include "physics/plasma/schemes/F10eViMdeltaf.h"
#include "physics/plasma/schemes/F10eF10iM.h"
#include "physics/plasma/schemes/F5eF10iM.h"
#include "physics/plasma/schemes/F5eF5iM.h"
#include "physics/plasma/schemes/MHD1TemperatureOhm.h"
#include "physics/plasma/criteria/criterion_position.h"
#include "physics/plasma/criteria/criterion_j.h"
#include "physics/plasma/criteria/criterion_j_ue.h"


namespace criterion {


int evaluate_criterion(model::Model* plasma_model[], model::Model* electromagnetic_model,
                       const Block_id& block_id, const Parameter& parameter) {

  int resulting_scheme_id = -1;

  if (parameter.criterion == nullptr) {
    resulting_scheme_id = parameter.default_scheme;
  } else if (parameter.criterion->get_criterion_id() == Parameter::kCriterionPosition) {
    resulting_scheme_id = ((criterion::Criterion_position*) parameter.criterion)->evaluate(block_id, parameter);
  } else if (parameter.criterion->get_criterion_id() == Parameter::kCriterionJ) {
    resulting_scheme_id = ((criterion::Criterion_j*) parameter.criterion)->evaluate(plasma_model, electromagnetic_model, parameter);
  } else if (parameter.criterion->get_criterion_id() == Parameter::kCriterionJUe) {
    resulting_scheme_id = ((criterion::Criterion_j_ue*) parameter.criterion)->evaluate(plasma_model, electromagnetic_model, parameter);
  }
  // else if (other criteria)
  else {
    std::cerr<<std::endl<<"Error: Criterion for scheme selection not found. "<<
      "Please check the value in parameter::criterion."<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool scheme_found_in_hierarchy = false;
  for (int s = 0; s < parameter.nschemes; ++s) {
    if (resulting_scheme_id == parameter.scheme_hierarchy[s]) {
      scheme_found_in_hierarchy = true;
      break;
    }
  }
  if (!scheme_found_in_hierarchy) {
    std::cerr<<std::endl<<"Error: The criterion requests a scheme that is not present "<<
      "in the parameter::scheme_hierarchy variable. Please adjust parameter::scheme_hierarchy "<<
      "or the criterion."<<std::endl<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  return resulting_scheme_id;
}


scheme::Scheme* allocate_scheme(int new_scheme_id, const Block_id& block_id,
        Mpi_boundary& boundary, const Parameter& parameter) {

  scheme::Scheme* selected_scheme = nullptr;

  switch (new_scheme_id) {
    case Parameter::kVeViM:
      selected_scheme = new scheme::VeViM(block_id, boundary, parameter);
      break;
    case Parameter::kF10eViM:
      selected_scheme = new scheme::F10eViM(block_id, boundary, parameter);
      break;
    case Parameter::kF10eF10iM:
      selected_scheme = new scheme::F10eF10iM(block_id, boundary, parameter);
      break;
    case Parameter::kF5eF10iM:
      selected_scheme = new scheme::F5eF10iM(block_id, boundary, parameter);
      break;
    case Parameter::kF5eF5iM:
      selected_scheme = new scheme::F5eF5iM(block_id, boundary, parameter);
      break;
    case Parameter::kMHD1TemperatureOhm:
      selected_scheme = new scheme::MHD1TemperatureOhm(block_id, boundary, parameter);
      break;
    case Parameter::kVeViP:
      selected_scheme = new scheme::VeViP(block_id, boundary, parameter);
      break;
    case Parameter::kVeViMdeltaf:
      selected_scheme = new scheme::VeViMdeltaf(block_id, boundary, parameter);
      break;
    case Parameter::kF10eViMdeltaf:
      selected_scheme = new scheme::F10eViMdeltaf(block_id, boundary, parameter);
      break;
    case Parameter::kVeViPdeltaf:
      selected_scheme = new scheme::VeViPdeltaf(block_id, boundary, parameter);
      break;
    default:
      std::cerr<<std::endl<<"Error: Selected scheme not found. "
        <<"Please check the value in parameter::default_scheme."<<std::endl;
  }

  return selected_scheme;
}


void convert_scheme(scheme::Scheme* current_scheme, scheme::Scheme* new_scheme, const Parameter& parameter) {

  model::Model* plasma_model[parameter.nspecies];
  model::Model* electromagnetic_model;

  current_scheme->calc_alloc_converted_model(electromagnetic_model, -1, new_scheme->get_scheme_id());
  for (int s = 0; s < parameter.nspecies; ++s) {
    current_scheme->calc_alloc_converted_model(plasma_model[s], s, new_scheme->get_scheme_id());
  }

  new_scheme->replace_dealloc_old_models(plasma_model, electromagnetic_model);
  
  new_scheme->set_dt(current_scheme->get_dt());

  if (parameter.electron_subcycling) {
    new_scheme->set_substeps(current_scheme->get_substeps());
    new_scheme->set_substep_counter(current_scheme->get_substep_counter());
  }
}


} // namespace criterion


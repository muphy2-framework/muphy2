/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "criterion_j.h"

#include <algorithm> // max_element
#include <math.h>

#include "framework/parameter.h"
#include "framework/utilities.h"
#include "framework/model.h"
#include "physics/plasma/models/vlasov.h"
#include "physics/plasma/models/fluid10.h"
#include "physics/plasma/models/fluid5.h"
#include "physics/plasma/models/mhd.h"
#include "physics/plasma/models/maxwell.h"
#include "physics/plasma/models/poisson.h"
#include "physics/plasma/models/ohm.h"
#include "physics/plasma/converters/converters.h"

namespace criterion {

Criterion_j::Criterion_j(int nused_schemes, int* scheme_hierarchy, real* thresholds) {

  nused_schemes_ = nused_schemes;

  scheme_hierarchy_ = new int[nused_schemes_];
  thresholds_ = new real[nused_schemes_-1];

  for (int i = 0; i < nused_schemes; ++i) {
    scheme_hierarchy_[i] = scheme_hierarchy[i];
  }
  for (int i = 0; i < nused_schemes-1; ++i) {
    thresholds_[i] = thresholds[i];
  }
}


Criterion_j::~Criterion_j() {
  delete[] thresholds_;
  delete[] scheme_hierarchy_;
}


int Criterion_j::get_criterion_id() { return Parameter::kCriterionJ; }


int Criterion_j::evaluate(model::Model* plasma_model[], model::Model* electromagnetic_model,
                           const Parameter& p) {

  real *j = new real[3 * p.ncells_x]();
  real *j_magnitude = new real[p.ncells_x]();
  real *un[p.nspecies];
  for (int i = 0; i < p.nspecies; ++i) {
    un[i] = new real[3 * p.ncells_x]();
  }
  #pragma acc enter data copyin(j[0:3*p.ncells_x],j_magnitude[0:p.ncells_x])
  #pragma acc enter data copyin(un[0:p.nspecies][0:3*p.ncells_x])

  for (int i = 0; i < p.nspecies; ++i) {

    switch (plasma_model[i]->get_model_id()) {
    case Parameter::kVlasov: {
      model::Vlasov * vlasov = (model::Vlasov *) plasma_model[i];
      convert::f_to_un(un[i], vlasov->f, p.species[i], p);
    } break;
    case Parameter::kFluid5: {
      model::Fluid5 * fluid5 = (model::Fluid5 *) plasma_model[i];
      util::copy_array(un[i], &fluid5->five_moments[p.ncells_x], 3*p.ncells_x);
    } break;
    case Parameter::kFluid10: {
      model::Fluid10 * fluid10 = (model::Fluid10 *) plasma_model[i];
      util::copy_array(un[i], &fluid10->ten_moments[p.ncells_x], 3*p.ncells_x);
    } break;
    case Parameter::kMHD1Temperature: {
      model::Ohm * ohm = (model::Ohm *)electromagnetic_model;
      j = ohm->J;
    } break;
      // other models
    }
  }

  if (plasma_model[0]->get_model_id() != Parameter::kMHD1Temperature) {
    convert::un_to_j(j, un[Parameter::kElectron], un[Parameter::kIon], p);
  }

  util::vector_magnitude(j_magnitude, j, p);

  #pragma acc update host(j_magnitude[0:p.ncells_x])
  real j_max = *std::max_element(j_magnitude,j_magnitude+p.ncells_x);

  #pragma acc exit data delete(un[0:p.nspecies][0:3*p.ncells_x])
  #pragma acc exit data delete(j[0:3*p.ncells_x],j_magnitude[0:p.ncells_x])
  for (int i = p.nspecies-1; i >= 0; --i) {
    delete[] un[i];
  }
  delete[] j_magnitude;
  delete[] j;


  for (int i = 0; i < nused_schemes_-1; ++i) {
    if (j_max >= thresholds_[i]) {
      return scheme_hierarchy_[i];
    }
  }
  return scheme_hierarchy_[nused_schemes_-1];

} // end evaluate()


} // namespace criterion


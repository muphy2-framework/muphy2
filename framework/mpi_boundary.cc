/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "mpi_boundary.h"
#include <vector>
#include <string.h>
#include "mpi.h"
#include "framework/definitions.h"
#include "framework/parameter.h"
#include "framework/utilities.h"


Mpi_boundary::Mpi_boundary(const MPI_Comm& comm3d, const int* neighbour_ranks,
            const int* my_coords, const Parameter& parameter)
: comm3d_(comm3d), neighbour_ranks_(neighbour_ranks), my_coords_(my_coords), p_(parameter)
{ }


Mpi_boundary::~Mpi_boundary() { }


real Mpi_boundary::mpi_min(real local_value)
{
  real global_min;
  MPI_Allreduce(&local_value, &global_min, 1, MPI_REAL_CXX, MPI_MIN, comm3d_);
  MPI_Barrier(comm3d_);

  return global_min;
}


int Mpi_boundary::mpi_min(int local_value)
{
  int global_min;
  MPI_Allreduce(&local_value, &global_min, 1, MPI_INT, MPI_MIN, comm3d_);
  MPI_Barrier(comm3d_);

  return global_min;
}

real Mpi_boundary::mpi_max(real local_value)
{
  real global_max;
  MPI_Allreduce(&local_value, &global_max, 1, MPI_REAL_CXX, MPI_MAX, comm3d_);
  MPI_Barrier(comm3d_);

  return global_max;
}


int Mpi_boundary::mpi_max(int local_value)
{
  int global_max;
  MPI_Allreduce(&local_value, &global_max, 1, MPI_INT, MPI_MAX, comm3d_);
  MPI_Barrier(comm3d_);

  return global_max;
}


real Mpi_boundary::mpi_sum(real local_value)
{
  real global_sum;
  MPI_Allreduce(&local_value, &global_sum, 1, MPI_REAL_CXX, MPI_SUM, comm3d_);
  MPI_Barrier(comm3d_);

  return global_sum;
}


int Mpi_boundary::mpi_sum(int local_value)
{
  int global_sum;
  MPI_Allreduce(&local_value, &global_sum, 1, MPI_INT, MPI_SUM, comm3d_);
  MPI_Barrier(comm3d_);

  return global_sum;
}


void Mpi_boundary::send_my_id(int my_id, int* neighbour_ids)
{
  for (int direction = 0; direction < p_.kNDirections; ++direction) {

    MPI_Irecv(&neighbour_ids[direction], 1, MPI_INT, neighbour_ranks_[direction], tags_[direction],
            comm3d_, &requests_[direction]);

    int tag_index = -1;
    switch (direction) {
      case 0: tag_index = 1; break;
      case 1: tag_index = 0; break;
      case 2: tag_index = 3; break;
      case 3: tag_index = 2; break;
      case 4: tag_index = 5; break;
      case 5: tag_index = 4; break;
    }

    MPI_Isend(&my_id, 1, MPI_INT, neighbour_ranks_[direction], tags_[tag_index],
            comm3d_, &requests_[direction+p_.kNDirections]);
  }

  MPI_Waitall(12, requests_, status_);
  MPI_Barrier(comm3d_);
}


void Mpi_boundary::send_neighbour_ids(int* neighbour_ids_send, int* neighbour_ids_receive)
{
  for (int direction = 0; direction < p_.kNDirections; ++direction) {

    MPI_Irecv(&neighbour_ids_receive[direction], 1, MPI_INT, neighbour_ranks_[direction], tags_[direction],
            comm3d_, &requests_[direction]);

    int tag_index = -1;
    switch (direction) {
      case 0: tag_index = 1; break;
      case 1: tag_index = 0; break;
      case 2: tag_index = 3; break;
      case 3: tag_index = 2; break;
      case 4: tag_index = 5; break;
      case 5: tag_index = 4; break;
    }

    MPI_Isend(&neighbour_ids_send[direction], 1, MPI_INT, neighbour_ranks_[direction], tags_[tag_index],
            comm3d_, &requests_[direction+p_.kNDirections]);
  }

  MPI_Waitall(12, requests_, status_);
  MPI_Barrier(comm3d_);
}


void Mpi_boundary::gather(int* gathered_array, int send)
{
  MPI_Gather(&send, 1, MPI_INT, gathered_array, 1, MPI_INT, 0, comm3d_);
}


void Mpi_boundary::gather(int* gathered_array, int* send, int send_array_length)
{
  MPI_Gather(send, send_array_length, MPI_INT, gathered_array, send_array_length, MPI_INT, 0, comm3d_);
}


void Mpi_boundary::scatter(int* scatter_array, int* receive)
{
  MPI_Scatter(scatter_array, 1, MPI_INT, receive, 1, MPI_INT, 0, comm3d_);
}


void Mpi_boundary::scatter(int* scatter_array, int* receive, int send_array_length)
{
  MPI_Scatter(scatter_array, send_array_length, MPI_INT, receive, send_array_length, MPI_INT, 0, comm3d_);
}


void Mpi_boundary::send_data(real* data, int array_size, int receiver_id)
{
  MPI_Isend(data, array_size, MPI_REAL_CXX, receiver_id, 0, comm3d_, &requests_[0]);
}


void Mpi_boundary::receive_data(real* data, int array_size, int sender_id)
{
  MPI_Irecv(data, array_size, MPI_REAL_CXX, sender_id, 0, comm3d_, &requests_[0]);
  MPI_Wait(&requests_[0], &status_[0]);
  MPI_Barrier(comm3d_);
}


int Mpi_boundary::calc_buffer_size(const int direction, const int ncomponents, const bool vspace, int species_id)
{
  int buffer_size = -1;

  if (direction == p_.kXL || direction == p_.kXU)
    buffer_size = ncomponents * p_.bd[0] * (p_.res_x[1] + 2*p_.bd[1]) * (p_.res_x[2] + 2*p_.bd[2]);
  else if (direction == p_.kYL || direction == p_.kYU)
    buffer_size = ncomponents * p_.bd[1] * (p_.res_x[0] + 2*p_.bd[0]) * (p_.res_x[2] + 2*p_.bd[2]);
  else if (direction == p_.kZL || direction == p_.kZU)
    buffer_size = ncomponents * p_.bd[2] * (p_.res_x[0] + 2*p_.bd[0]) * (p_.res_x[1] + 2*p_.bd[1]);

  if (vspace && buffer_size != -1)
    buffer_size *= p_.species[species_id].res_v[0]*p_.species[species_id].res_v[1]*p_.species[species_id].res_v[2];

  return buffer_size;
}


int* Mpi_boundary::calc_buffer_array_dimensions(const int direction)
{
  static int buffer_array_dimensions[3];

  if (direction == p_.kXL || direction == p_.kXU) {
    buffer_array_dimensions[0] = p_.bd[0];
    buffer_array_dimensions[1] = p_.res_x[1] + 2*p_.bd[1];
    buffer_array_dimensions[2] = p_.res_x[2] + 2*p_.bd[2];
  }
  else if (direction == p_.kYL || direction == p_.kYU) {
    buffer_array_dimensions[0] = p_.res_x[0] + 2*p_.bd[0];
    buffer_array_dimensions[1] = p_.bd[1];
    buffer_array_dimensions[2] = p_.res_x[2] + 2*p_.bd[2];
  }
  else if (direction == p_.kZL || direction == p_.kZU) {
    buffer_array_dimensions[0] = p_.res_x[0] + 2*p_.bd[0];
    buffer_array_dimensions[1] = p_.res_x[1] + 2*p_.bd[1];
    buffer_array_dimensions[2] = p_.bd[2];
  }

  return buffer_array_dimensions;
}


void Mpi_boundary::exchange_buffer(real* buffer_out, real* buffer_in, const int buffer_size, const int direction)
{
  if ((p_.dimensionality_x == 2 && direction >= 4) || (p_.dimensionality_x == 1 && direction >= 2))
    return;

  #pragma acc host_data use_device(buffer_in)
  MPI_Irecv(buffer_in, buffer_size, MPI_REAL_CXX, neighbour_ranks_[direction], tags_[direction],
          comm3d_, &requests_[direction]);

  int tag_index = -1;
  switch (direction) {
    case 0: tag_index = 1; break;
    case 1: tag_index = 0; break;
    case 2: tag_index = 3; break;
    case 3: tag_index = 2; break;
    case 4: tag_index = 5; break;
    case 5: tag_index = 4; break;
  }

  #pragma acc host_data use_device(buffer_out)
  MPI_Isend(buffer_out, buffer_size, MPI_REAL_CXX, neighbour_ranks_[direction], tags_[tag_index],
          comm3d_, &requests_[direction+p_.kNDirections]);
}


void Mpi_boundary::wait(const bool* wait_for_this_direction)
{
  std::vector<MPI_Request> requests;  

  for (int direction = 0; direction < p_.kNDirections; ++direction) {

    if ((p_.dimensionality_x == 2 && direction >= 4) ||
        (p_.dimensionality_x == 1 && direction >= 2)) {
      break;
    }

    if (wait_for_this_direction[direction]) {
      requests.push_back(requests_[direction]);
      requests.push_back(requests_[direction+p_.kNDirections]);
    }
  }
  if (requests.size() > 0) {
    MPI_Status status[requests.size()];
    MPI_Waitall(requests.size(), requests.data(), status);
  }
}


void Mpi_boundary::wait()
{
  bool wait_for_this_direction[6] = {1, 1, 1, 1, 1, 1};
  wait(wait_for_this_direction);
  // barrier for robustness
  // (possible because all directions are exchanged)
  MPI_Barrier(comm3d_);
}


void Mpi_boundary::extract_buffer(const real* u, real* buffer,
        const int direction, const int ncomponents, const bool vspace, int species_id)
{
  if ((p_.dimensionality_x == 2 && direction >= 4) || (p_.dimensionality_x == 1 && direction >= 2))
    return;

  if (!vspace) {
    switch (direction) {
      case Parameter::kXL: extract_buffer_x(u, buffer, ncomponents, 0); break;
      case Parameter::kXU: extract_buffer_x(u, buffer, ncomponents, 1); break;
      case Parameter::kYL: extract_buffer_y(u, buffer, ncomponents, 0); break;
      case Parameter::kYU: extract_buffer_y(u, buffer, ncomponents, 1); break;
      case Parameter::kZL: extract_buffer_z(u, buffer, ncomponents, 0); break;
      case Parameter::kZU: extract_buffer_z(u, buffer, ncomponents, 1); break;
    }
  }
  else {
    switch (direction) {
      case Parameter::kXL: extract_buffer_f_x(u, buffer, 0, species_id); break;
      case Parameter::kXU: extract_buffer_f_x(u, buffer, 1, species_id); break;
      case Parameter::kYL: extract_buffer_f_y(u, buffer, 0, species_id); break;
      case Parameter::kYU: extract_buffer_f_y(u, buffer, 1, species_id); break;
      case Parameter::kZL: extract_buffer_f_z(u, buffer, 0, species_id); break;
      case Parameter::kZU: extract_buffer_f_z(u, buffer, 1, species_id); break;
    }
  }
}


void Mpi_boundary::apply_buffer(real* u, const real* buffer,
        const int direction, const int ncomponents, const bool vspace, int species_id)
{
  if (!vspace) {
    if ((p_.dimensionality_x == 2 && direction >= 4) || (p_.dimensionality_x == 1 && direction >= 2)) {
      dummy_dimension_handling(u, ncomponents, direction);
    } else {
      switch (direction) {
        case Parameter::kXL: apply_buffer_x(u, buffer, ncomponents, 0); break;
        case Parameter::kXU: apply_buffer_x(u, buffer, ncomponents, 1); break;
        case Parameter::kYL: apply_buffer_y(u, buffer, ncomponents, 0); break;
        case Parameter::kYU: apply_buffer_y(u, buffer, ncomponents, 1); break;
        case Parameter::kZL: apply_buffer_z(u, buffer, ncomponents, 0); break;
        case Parameter::kZU: apply_buffer_z(u, buffer, ncomponents, 1); break;
      }
    }
  }
  else {
    if ((p_.dimensionality_x == 2 && direction >= 4) || (p_.dimensionality_x == 1 && direction >= 2))
      dummy_dimension_handling_f(u, direction, species_id);
    else {
      switch (direction) {
        case Parameter::kXL: apply_buffer_f_x(u, buffer, 0, species_id); break;
        case Parameter::kXU: apply_buffer_f_x(u, buffer, 1, species_id); break;
        case Parameter::kYL: apply_buffer_f_y(u, buffer, 0, species_id); break;
        case Parameter::kYU: apply_buffer_f_y(u, buffer, 1, species_id); break;
        case Parameter::kZL: apply_buffer_f_z(u, buffer, 0, species_id); break;
        case Parameter::kZU: apply_buffer_f_z(u, buffer, 1, species_id); break;
      }
    }
  }
}


void Mpi_boundary::apply_buffer_schwarz(real* u, const real* buffer, const int direction)
{
  if ((p_.dimensionality_x == 2 && direction >= 4) || (p_.dimensionality_x == 1 && direction >= 2))
    dummy_dimension_handling(u, 1, direction);
  else {
    switch (direction) {
      case Parameter::kXL: apply_buffer_schwarz_x(u, buffer, 0); break;
      case Parameter::kXU: apply_buffer_schwarz_x(u, buffer, 1); break;
      case Parameter::kYL: apply_buffer_schwarz_y(u, buffer, 0); break;
      case Parameter::kYU: apply_buffer_schwarz_y(u, buffer, 1); break;
      case Parameter::kZL: apply_buffer_schwarz_z(u, buffer, 0); break;
      case Parameter::kZU: apply_buffer_schwarz_z(u, buffer, 1); break;
    }
  }
}


void Mpi_boundary::apply_boundary_conditions(real* u, const std::string* bd_cond, const std::string* orientation,
        const int direction, const int ncomponents, const bool vspace, int species_id)
{
  // return if block is not at domain border
  if ((direction == Parameter::kXL && my_coords_[0] != 0) ||
      (direction == Parameter::kXU && my_coords_[0] != p_.nproc[0]-1) ||
      (direction == Parameter::kYL && my_coords_[1] != 0) ||
      (direction == Parameter::kYU && my_coords_[1] != p_.nproc[1]-1) ||
      (direction == Parameter::kZL && my_coords_[2] != 0) ||
      (direction == Parameter::kZU && my_coords_[2] != p_.nproc[2]-1)) {
    return;
  }

  if (!vspace) {
    switch (direction) {
      case Parameter::kXL: apply_boundary_conditions_x(u, ncomponents, 0, bd_cond, orientation); break;
      case Parameter::kXU: apply_boundary_conditions_x(u, ncomponents, 1, bd_cond, orientation); break;
      case Parameter::kYL: apply_boundary_conditions_y(u, ncomponents, 0, bd_cond, orientation); break;
      case Parameter::kYU: apply_boundary_conditions_y(u, ncomponents, 1, bd_cond, orientation); break;
      case Parameter::kZL: apply_boundary_conditions_z(u, ncomponents, 0, bd_cond, orientation); break;
      case Parameter::kZU: apply_boundary_conditions_z(u, ncomponents, 1, bd_cond, orientation); break;
    }
  }
  else {
    switch (direction) {
      case Parameter::kXL: apply_boundary_conditions_f_x(u, 0, bd_cond[0], species_id); break;
      case Parameter::kXU: apply_boundary_conditions_f_x(u, 1, bd_cond[0], species_id); break;
      case Parameter::kYL: apply_boundary_conditions_f_y(u, 0, bd_cond[1], species_id); break;
      case Parameter::kYU: apply_boundary_conditions_f_y(u, 1, bd_cond[1], species_id); break;
      case Parameter::kZL: apply_boundary_conditions_f_z(u, 0, bd_cond[2], species_id); break;
      case Parameter::kZU: apply_boundary_conditions_f_z(u, 1, bd_cond[2], species_id); break;
    }
  }
}


void Mpi_boundary::exchange_field(real* u, const std::string* bd_cond, const std::string* orientation,
          const bool* exchange_in_this_direction, const int ncomponents)
{ // ncomponents = 1 for scalar, = 3 for vector, etc.

  real* buffer_out[Parameter::kNDirections];
  real* buffer_in [Parameter::kNDirections];
  int buffer_size[Parameter::kNDirections];

  for (int d = 0; d < Parameter::kNDirections; ++d) {
    if (exchange_in_this_direction[d]) {
      buffer_size[d] = calc_buffer_size(d, ncomponents);
      buffer_out[d] = new real[buffer_size[d]];
      buffer_in[d]  = new real[buffer_size[d]];
      #pragma acc enter data create(buffer_out[d:1][0:buffer_size[d]],\
                                    buffer_in [d:1][0:buffer_size[d]])

      extract_buffer(u, buffer_out[d], d, ncomponents);
      exchange_buffer(buffer_out[d], buffer_in[d], buffer_size[d], d);
    }
  }

  wait(exchange_in_this_direction);

  for (int d = 0; d < Parameter::kNDirections; ++d) {
    if (exchange_in_this_direction[d]) {
      apply_buffer(u, buffer_in[d], d, ncomponents);

      #pragma acc exit data delete(buffer_out[d:1][0:buffer_size[d]],\
                                   buffer_in [d:1][0:buffer_size[d]])
      delete[] buffer_out[d];
      delete[] buffer_in [d];
    }
    apply_boundary_conditions(u, bd_cond, orientation, d, ncomponents);
  }
}


void Mpi_boundary::exchange_field(real* u, const std::string* bd_cond, const std::string* orientation,
        const int ncomponents)
{
  bool exchange_in_this_direction[6] = {1, 1, 1, 1, 1, 1};
  exchange_field(u, bd_cond, orientation, exchange_in_this_direction, ncomponents);
}


void Mpi_boundary::exchange_schwarz(real* u, const std::string* bd_cond, const std::string* orientation)
{
  real* buffer_out[Parameter::kNDirections];
  real* buffer_in [Parameter::kNDirections];
  int buffer_size[Parameter::kNDirections];

  for (int d = 0; d < Parameter::kNDirections; ++d) {
    buffer_size[d] = calc_buffer_size(d, 1);
    buffer_out[d] = new real[buffer_size[d]];
    buffer_in[d]  = new real[buffer_size[d]];
    #pragma acc enter data create(buffer_out[d:1][0:buffer_size[d]],\
                                  buffer_in [d:1][0:buffer_size[d]])

    extract_buffer(u, buffer_out[d], d, 1);
    exchange_buffer(buffer_out[d], buffer_in[d], buffer_size[d], d);
  }

  wait();

  for (int d = 0; d < Parameter::kNDirections; ++d) {
    apply_buffer_schwarz(u, buffer_in[d], d);
    apply_boundary_conditions(u, bd_cond, orientation, d, 1);

    #pragma acc exit data delete(buffer_out[d:1][0:buffer_size[d]],\
                                 buffer_in [d:1][0:buffer_size[d]])
    delete[] buffer_out[d];
    delete[] buffer_in [d];
  }
}


void Mpi_boundary::exchange_distribution_function(real* f, const std::string* bd_cond,
        const bool* exchange_in_this_direction)
{
  real* buffer_out[Parameter::kNDirections];
  real* buffer_in [Parameter::kNDirections];
  int buffer_size[Parameter::kNDirections];

  for (int d = 0; d < Parameter::kNDirections; ++d) {

    if (exchange_in_this_direction[d]) {
      buffer_size[d] = calc_buffer_size(d, 1, true);
      buffer_out[d] = new real[buffer_size[d]];
      buffer_in[d]  = new real[buffer_size[d]];
      #pragma acc enter data create(buffer_out[d:1][0:buffer_size[d]],\
                                    buffer_in [d:1][0:buffer_size[d]])

      extract_buffer(f, buffer_out[d], d, 1, true);
      exchange_buffer(buffer_out[d], buffer_in[d], buffer_size[d], d);
    }
  }

  wait(exchange_in_this_direction);

  for (int d = 0; d < Parameter::kNDirections; ++d) {
    std::string orientation[1] = {"ccc"};
    if (exchange_in_this_direction[d]) {
      apply_buffer(f, buffer_in[d], d, 1, true);

      #pragma acc exit data delete(buffer_out[d:1][0:buffer_size[d]],\
                                   buffer_in [d:1][0:buffer_size[d]])
      delete[] buffer_out[d];
      delete[] buffer_in [d];
    }
    apply_boundary_conditions(f, bd_cond, orientation, d, 1, true);
  }
}


void Mpi_boundary::exchange_distribution_function(real* f, const std::string* bd_cond)
{
  bool exchange_in_this_direction[6] = {1, 1, 1, 1, 1, 1};
  exchange_distribution_function(f, bd_cond, exchange_in_this_direction);
}


// fortran wrappers
void Mpi_boundary::dummy_dimension_handling(real* u, const int ncomponents, const int direction)
{
  boundary_dummy_dimension_handling_(u, direction, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::extract_buffer_x(const real* u, real* out_buffer, const int ncomponents, const int side)
{
  boundary_extract_buffer_x_(u, out_buffer, side, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::extract_buffer_y(const real* u, real* out_buffer, const int ncomponents, const int side)
{
  boundary_extract_buffer_y_(u, out_buffer, side, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::extract_buffer_z(const real* u, real* out_buffer, const int ncomponents, const int side)
{
  boundary_extract_buffer_z_(u, out_buffer, side, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::apply_buffer_x(real* u, const real* in_buffer, const int ncomponents, const int side)
{
  boundary_apply_buffer_x_(u, in_buffer, side, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::apply_buffer_y(real* u, const real* in_buffer, const int ncomponents, const int side)
{
  boundary_apply_buffer_y_(u, in_buffer, side, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::apply_buffer_z(real* u, const real* in_buffer, const int ncomponents, const int side)
{
  boundary_apply_buffer_z_(u, in_buffer, side, p_.res_x_minus_one, p_.bd, ncomponents);
}

void Mpi_boundary::apply_buffer_schwarz_x(real* u, const real* in_buffer, const int side)
{
  boundary_apply_buffer_schwarz_x_(u, in_buffer, side, p_.res_x_minus_one, p_.bd);
}

void Mpi_boundary::apply_buffer_schwarz_y(real* u, const real* in_buffer, const int side)
{
  boundary_apply_buffer_schwarz_y_(u, in_buffer, side, p_.res_x_minus_one, p_.bd);
}

void Mpi_boundary::apply_buffer_schwarz_z(real* u, const real* in_buffer, const int side)
{
  boundary_apply_buffer_schwarz_z_(u, in_buffer, side, p_.res_x_minus_one, p_.bd);
}

void Mpi_boundary::apply_boundary_conditions_x(real* u, const int ncomponents, const int side,
        const std::string bd_cond[], const std::string orientation[])
{
  char bd_cond_fortran[6*ncomponents], orientation_fortran[3*ncomponents];
  for (int i = 0; i < ncomponents; ++i) {
    util::str_cpp_to_fortran(bd_cond[i], &bd_cond_fortran[i*6]);
    util::str_cpp_to_fortran(orientation[i], &orientation_fortran[i*3]);
  }

  boundary_apply_boundary_conditions_x_(u, side, p_.res_x_minus_one,
          p_.bd, ncomponents, bd_cond_fortran, orientation_fortran);
}

void Mpi_boundary::apply_boundary_conditions_y(real* u, const int ncomponents, const int side,
        const std::string bd_cond[], const std::string orientation[])
{
  char bd_cond_fortran[6*ncomponents], orientation_fortran[3*ncomponents];
  for (int i = 0; i < ncomponents; ++i) {
    util::str_cpp_to_fortran(bd_cond[i], &bd_cond_fortran[i*6]);
    util::str_cpp_to_fortran(orientation[i], &orientation_fortran[i*3]);
  }

  boundary_apply_boundary_conditions_y_(u, side, p_.res_x_minus_one,
          p_.bd, ncomponents, bd_cond_fortran, orientation_fortran);
}

void Mpi_boundary::apply_boundary_conditions_z(real* u, const int ncomponents, const int side,
        const std::string bd_cond[], const std::string orientation[])
{
  char bd_cond_fortran[6*ncomponents], orientation_fortran[3*ncomponents];
  for (int i = 0; i < ncomponents; ++i) {
    util::str_cpp_to_fortran(bd_cond[i], &bd_cond_fortran[i*6]);
    util::str_cpp_to_fortran(orientation[i], &orientation_fortran[i*3]);
  }

  boundary_apply_boundary_conditions_z_(u, side, p_.res_x_minus_one,
          p_.bd, ncomponents, bd_cond_fortran, orientation_fortran);
}

void Mpi_boundary::dummy_dimension_handling_f(real* f, const int direction, int species_id)
{
  boundary_dummy_dimension_handling_f_(f, direction, p_.res_x_minus_one,
                                       p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::extract_buffer_f_x(const real* f, real* out_buffer, const int side, int species_id)
{
  boundary_extract_buffer_f_x_(f, out_buffer, side, p_.res_x_minus_one, p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::extract_buffer_f_y(const real* f, real* out_buffer, const int side, int species_id)
{
  boundary_extract_buffer_f_y_(f, out_buffer, side, p_.res_x_minus_one, p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::extract_buffer_f_z(const real* f, real* out_buffer, const int side, int species_id)
{
  boundary_extract_buffer_f_z_(f, out_buffer, side, p_.res_x_minus_one, p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::apply_buffer_f_x(real* f, const real* in_buffer, const int side, int species_id)
{
  boundary_apply_buffer_f_x_(f, in_buffer, side, p_.res_x_minus_one, p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::apply_buffer_f_y(real* f, const real* in_buffer, const int side, int species_id)
{
  boundary_apply_buffer_f_y_(f, in_buffer, side, p_.res_x_minus_one, p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::apply_buffer_f_z(real* f, const real* in_buffer, const int side, int species_id)
{
  boundary_apply_buffer_f_z_(f, in_buffer, side, p_.res_x_minus_one, p_.species[species_id].res_v, p_.bd);
}

void Mpi_boundary::apply_boundary_conditions_f_x(real* f, const int side,
        const std::string bd_cond, int species_id)
{
  char bd_cond_fortran[6];
  util::str_cpp_to_fortran(bd_cond, bd_cond_fortran);

  boundary_apply_boundary_conditions_f_x_(f, side, p_.res_x_minus_one,
          p_.species[species_id].res_v, p_.bd, bd_cond_fortran);
}

void Mpi_boundary::apply_boundary_conditions_f_y(real* f, const int side,
        const std::string bd_cond, int species_id)
{
  char bd_cond_fortran[6];
  util::str_cpp_to_fortran(bd_cond, bd_cond_fortran);

  boundary_apply_boundary_conditions_f_y_(f, side, p_.res_x_minus_one,
          p_.species[species_id].res_v, p_.bd, bd_cond_fortran);
}

void Mpi_boundary::apply_boundary_conditions_f_z(real* f, const int side,
        const std::string bd_cond, int species_id)
{
  char bd_cond_fortran[6];
  util::str_cpp_to_fortran(bd_cond, bd_cond_fortran);

  boundary_apply_boundary_conditions_f_z_(f, side, p_.res_x_minus_one,
          p_.species[species_id].res_v, p_.bd, bd_cond_fortran);
}


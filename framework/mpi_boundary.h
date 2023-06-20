/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <sstream>
#include <string.h>
#include "mpi.h"
#include "framework/definitions.h"

struct Parameter;

// fortran functions
extern "C" void boundary_dummy_dimension_handling_(real* u, const int& direction,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_extract_buffer_x_(const real* u, real* out_buffer, const int& side,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_extract_buffer_y_(const real* u, real* out_buffer, const int& side,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_extract_buffer_z_(const real* u, real* out_buffer, const int& side,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_apply_buffer_x_(real* u, const real* in_buffer, const int& side,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_apply_buffer_y_(real* u, const real* in_buffer, const int& side,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_apply_buffer_z_(real* u, const real* in_buffer, const int& side,
          const int* dimX, const int* BD, const int& ncomponents);
extern "C" void boundary_apply_buffer_schwarz_x_(real* u, const real* in_buffer, const int& side,
          const int* dimX, const int* BD);
extern "C" void boundary_apply_buffer_schwarz_y_(real* u, const real* in_buffer, const int& side,
          const int* dimX, const int* BD);
extern "C" void boundary_apply_buffer_schwarz_z_(real* u, const real* in_buffer, const int& side,
          const int* dimX, const int* BD);
extern "C" void boundary_apply_boundary_conditions_x_(real* u, const int& side,
          const int* dimX, const int* BD, const int& ncomponents, const char* bdCond, const char* orientation);
extern "C" void boundary_apply_boundary_conditions_y_(real* u, const int& side,
          const int* dimX, const int* BD, const int& ncomponents, const char* bdCond, const char* orientation);
extern "C" void boundary_apply_boundary_conditions_z_(real* u, const int& side,
          const int* dimX, const int* BD, const int& ncomponents, const char* bdCond, const char* orientation);

extern "C" void boundary_dummy_dimension_handling_f_(real* f, const int& direction,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_extract_buffer_f_x_(const real* f, real* out_buffer, const int& side,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_extract_buffer_f_y_(const real* f, real* out_buffer, const int& side,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_extract_buffer_f_z_(const real* f, real* out_buffer, const int& side,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_apply_buffer_f_x_(real* f, const real* in_buffer, const int& side,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_apply_buffer_f_y_(real* f, const real* in_buffer, const int& side,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_apply_buffer_f_z_(real* f, const real* in_buffer, const int& side,
          const int* dimX, const int* dimV, const int* BD);
extern "C" void boundary_apply_boundary_conditions_f_x_(real* f, const int& side,
          const int* dimX, const int* dimV, const int* BD, const char* bdCond);
extern "C" void boundary_apply_boundary_conditions_f_y_(real* f, const int& side,
          const int* dimX, const int* dimV, const int* BD, const char* bdCond);
extern "C" void boundary_apply_boundary_conditions_f_z_(real* f, const int& side,
          const int* dimX, const int* dimV, const int* BD, const char* bdCond);

class Mpi_boundary {

  public:
    Mpi_boundary(const MPI_Comm& comm3d, const int* neighbour_ranks,
            const int* my_coords, const Parameter& parameter);
    ~Mpi_boundary();

    real mpi_min(real local_value);
    int mpi_min(int local_value);
    real mpi_max(real local_value);
    int mpi_max(int local_value);
    real mpi_sum(real local_value);
    int mpi_sum(int local_value);

    // wrappers for fortran functions
    void dummy_dimension_handling(real* u, const int ncomponents, const int direction);
    void extract_buffer_x(const real* u, real* out_buffer, const int ncomponents, const int side);
    void extract_buffer_y(const real* u, real* out_buffer, const int ncomponents, const int side);
    void extract_buffer_z(const real* u, real* out_buffer, const int ncomponents, const int side);
    void apply_buffer_x(real* u, const real* in_buffer, const int ncomponents, const int side);
    void apply_buffer_y(real* u, const real* in_buffer, const int ncomponents, const int side);
    void apply_buffer_z(real* u, const real* in_buffer, const int ncomponents, const int side);
    void apply_buffer_schwarz_x(real* u, const real* in_buffer, const int side);
    void apply_buffer_schwarz_y(real* u, const real* in_buffer, const int side);
    void apply_buffer_schwarz_z(real* u, const real* in_buffer, const int side);
    void apply_boundary_conditions_x(real* u, const int ncomponents, const int side,
             const std::string bd_cond[], const std::string orientation[]);
    void apply_boundary_conditions_y(real* u, const int ncomponents, const int side,
             const std::string bd_cond[], const std::string orientation[]);
    void apply_boundary_conditions_z(real* u, const int ncomponents, const int side,
             const std::string bd_cond[], const std::string orientation[]);

    void dummy_dimension_handling_f(real* f, int direction, int species_id);
    void extract_buffer_f_x(const real* f, real* out_buffer, const int side, int species_id);
    void extract_buffer_f_y(const real* f, real* out_buffer, const int side, int species_id);
    void extract_buffer_f_z(const real* f, real* out_buffer, const int side, int species_id);
    void apply_buffer_f_x(real* f, const real* in_buffer, const int side, int species_id);
    void apply_buffer_f_y(real* f, const real* in_buffer, const int side, int species_id);
    void apply_buffer_f_z(real* f, const real* in_buffer, const int side, int species_id);
    void apply_boundary_conditions_f_x(real* f, const int side, const std::string bd_cond, int species_id);
    void apply_boundary_conditions_f_y(real* f, const int side, const std::string bd_cond, int species_id);
    void apply_boundary_conditions_f_z(real* f, const int side, const std::string bd_cond, int species_id);

    // boundary routines
    void exchange_buffer(real* buffer_out, real* buffer_in, const int buffer_size, const int direction);
    void wait();
    void wait(const bool* wait_for_this_direction);
    int calc_buffer_size(const int direction, const int ncomponents, const bool vspace=false, int species_id=0);
    int* calc_buffer_array_dimensions(const int direction);
    void extract_buffer(const real* u, real* buffer,
            const int direction, const int ncomponents, const bool vspace=false, int species_id=0);
    void apply_buffer(real* u, const real* buffer,
            const int direction, const int ncomponents, const bool vspace=false, int species_id=0);
    void apply_buffer_schwarz(real* u, const real* buffer, const int direction);
    void apply_boundary_conditions(real* u, const std::string* bd_cond, const std::string* orientation,
            const int direction, const int ncomponents, const bool vspace=false, int species_id=0);
    void send_my_id(int my_id, int* neighbour_ids);
    void send_neighbour_ids(int* neighbour_ids_send, int* neighbour_ids_receive);
    void gather(int* gathered_array, int send);
    void gather(int* gathered_array, int* send, int send_array_length);
    void scatter(int* scatter_array, int* receive);
    void scatter(int* scatter_array, int* receive, int send_array_length);
    void send_data(real* data, int array_size, int receiver_id);
    void receive_data(real* data, int array_size, int sender_id);
    void exchange_field(real* u, const std::string* bd_cond,
            const std::string* orientation, const bool* exchange_in_this_direction,
            const int ncomponents);
    void exchange_field(real* u, const std::string* bd_cond,
            const std::string* orientation, const int ncomponents);
    void exchange_schwarz(real* u, const std::string* bd_cond, const std::string* orientation);
    void exchange_distribution_function(real* f, const std::string* bd_cond,
            const bool* exchange_in_this_direction);
    void exchange_distribution_function(real* f, const std::string* bd_cond);

    //int max_local_rank;

  private:
    const MPI_Comm& comm3d_;
    const int* neighbour_ranks_,
             * my_coords_;

    const Parameter& p_;

    int tags_[6] = {0,1,2,3,4,5};
    MPI_Status status_[12];
    MPI_Request requests_[12];
};

#endif // BOUNDARY_H_

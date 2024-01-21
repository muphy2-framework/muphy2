/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "block.h"

#include "mpi.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm> // std::min_element
#include <string.h>
#include <sys/stat.h> // mkdir etc
#include <cstdlib> // exit
#include <chrono>

#include "framework/definitions.h"
#include "framework/restart.h"
#include "framework/mpi_boundary.h"
#include "framework/parameter.h"
#include "framework/output_write.h"
#include "framework/scheme.h"
#include "physics/plasma/criteria/select_scheme.h"

Block::Block()
: boundary_(comm3d_, block_id_.neighbour_ids, block_id_.my_coords, parameter_)
{
  // MPI initialization
  int procs;
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &block_id_.my_id);

  MPI_Comm node_comm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
  int local_rank = -1;
  MPI_Comm_rank(node_comm, &local_rank);
  block_id_.selected_device = local_rank;
  if (parameter_.nslots_per_node > 0) {
    block_id_.selected_device = local_rank%parameter_.nslots_per_node;
  }
  #pragma acc set device_num(block_id_.selected_device)
  MPI_Comm_free(&node_comm);
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int processor_name_string_length;
  MPI_Get_processor_name(processor_name, &processor_name_string_length);

  if (parameter_.nproc[0]*parameter_.nproc[1]*parameter_.nproc[2] != procs) {
    std::cerr  << std::endl << "Error: Number of MPI processes not compatible with MPI_Cart_create layout " <<
      "defined in parameter::nproc." << std::endl <<
      "Number of processes launched by mpirun = " << procs << "," << std::endl <<
      "nproc[0]*nproc[1]*nproc[2] = " << parameter_.nproc[0]*parameter_.nproc[1]*parameter_.nproc[2] <<
      std::endl << "(must be the same)." << std::endl << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int periods[3] = {1,1,1};
  int reorder = 0;
  MPI_Cart_create(MPI_COMM_WORLD, 3, parameter_.nproc, periods, reorder, &comm3d_);
  MPI_Cart_coords(comm3d_, block_id_.my_id, 3, block_id_.my_coords);
  MPI_Cart_shift(comm3d_, 0, 1, &block_id_.neighbour_ids[0], &block_id_.neighbour_ids[1]);
  MPI_Cart_shift(comm3d_, 1, 1, &block_id_.neighbour_ids[2], &block_id_.neighbour_ids[3]);
  MPI_Cart_shift(comm3d_, 2, 1, &block_id_.neighbour_ids[4], &block_id_.neighbour_ids[5]);
  //

  max_local_rank_ = boundary_.mpi_max(local_rank);

  const int nprocs = parameter_.nproc[0]*parameter_.nproc[1]*parameter_.nproc[2],
            nnodes = nprocs/(max_local_rank_+1);
  block_id_.compute_node = floor(nnodes * ((real)block_id_.my_id)/nprocs);

  if (block_id_.my_id == 0) {
    std::cout << std::endl << "muphy 2.0"<< std::endl
    << "Multiphysics plasma simulation framework"<< std::endl
    << std::endl
    << "----------------------------------------------------------------------------" << std::endl
    << " If this software contributes to findings you decide to present or publish,"  << std::endl
    << " please be so kind and cite the reference below. Thank you!"                  << std::endl
    << std::endl
    << " Allmann-Rahn, F., Lautenbach, S., Deisenhofer, M., Grauer, R., 2024."        << std::endl 
    << " The muphyII code: Multiphysics plasma simulation on large HPC systems."      << std::endl 
    << " Computer Physics Communications. https://doi.org/10.1016/j.cpc.2023.109064"  << std::endl
    << "----------------------------------------------------------------------------" << std::endl;

    // create output directory
    int dir_err = mkdir(parameter_.output_directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (dir_err < 0 && errno != EEXIST) {
      std::cerr<<std::endl<<"Error: Could not create output directory. "<<
        "Only the highest level directory is created. Please make sure that "<<
        "the parent directory exists and/or check write permissions."<<
        std::endl<<std::endl;
      std::exit(EXIT_FAILURE);
    }
    mkdir((parameter_.output_directory + "/vti").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (parameter_.phase_space_output != "none") {
      mkdir((parameter_.output_directory + "/phase_space").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    if (parameter_.phase_space_output == "slice") {
      mkdir((parameter_.output_directory + "/phase_space/vti").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    mkdir((parameter_.output_directory + "/csv").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((parameter_.output_directory + "/log").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (parameter_.restart_output) {
      mkdir((parameter_.output_directory + "/restart").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    // TODO replace system() calls by safer alternative
    if (!parameter_.restart) {
      system(("rm " + parameter_.output_directory + "/csv/*.csv >>/dev/null 2>>/dev/null").c_str());
    }

    parameter_.write_initial_conditions_files();
    
    std::cout << std::endl << "Output directory: " << parameter_.output_directory
              << std::endl << std::endl << "Initializing ..." << std::endl;
  }

  for (int i = 0; i < 3; ++i) {
    parameter_.xb_loc[i] = parameter_.xb[i] + block_id_.my_coords[i] * parameter_.res_x[i] * parameter_.dx[i];
  }

  dt_update_scheme_ = (parameter_.t_end - parameter_.t_begin) / parameter_.nupdates_scheme;
  dt_output_vtk_ = (parameter_.t_end - parameter_.t_begin) / parameter_.noutputs_vtk;
  dt_output_csv_ = (parameter_.t_end - parameter_.t_begin) / parameter_.noutputs_csv;

  freeze_dt_default_ = parameter_.freeze_dt; // for disabling dt change during consecutive outputs

  if (parameter_.restart) {
    output_number_ = get_t_begin() / dt_output_vtk_ + 1;
  }

  // set up scheme
  int selected_scheme = -1;
  
  if (parameter_.criterion != nullptr) {
    // check if default scheme is present in scheme hierarchy
    bool default_scheme_found_in_hierarchy = false;
    for (int s = 0; s < parameter_.nschemes; ++s) {
      if (parameter_.default_scheme == parameter_.scheme_hierarchy[s]) {
        default_scheme_found_in_hierarchy = true;
        break;
      }
    }
    if (!default_scheme_found_in_hierarchy) {
      std::cerr<<std::endl<<"Error: The default scheme must be present in the "<<
        "scheme hierarchy. Please adust either parameter::default_scheme "<<
        "or parameter::scheme_hierarchy."<<std::endl<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  if (parameter_.restart) {
    selected_scheme = restart::load_scheme_id(parameter_.output_directory, block_id_.my_coords);
  } else if (parameter_.criterion == nullptr) {
    selected_scheme = parameter_.default_scheme;
  } else {
    scheme::Scheme *scheme_tmp = criterion::allocate_scheme(parameter_.default_scheme,
                                 block_id_, boundary_, parameter_);
    scheme_tmp->init();
    selected_scheme = scheme_tmp->evaluate_criterion();
    delete scheme_tmp;
  }
  scheme_ = criterion::allocate_scheme(selected_scheme, block_id_, boundary_, parameter_);
  block_id_.my_scheme_id = scheme_->get_scheme_id();
  boundary_.send_my_id(block_id_.my_scheme_id, block_id_.neighbour_scheme_ids);
  scheme_->init();

  if (parameter_.nslots_per_node > 0 && parameter_.criterion != nullptr) {
    load_balancing();
  }

  // initial output
  if (!parameter_.restart) {
    if (parameter_.noutputs_vtk != 0) {
      scheme_->output(parameter_.t_end, output_number_, true);
      ++output_number_;
    } else if (parameter_.noutputs_csv != 0) {
      scheme_->output(parameter_.t_end, output_number_, false);
    }
  }
}


Block::~Block() {

  if (block_id_.my_id == 0) {
    std::cout << "Run completed." << std::endl << std::endl;
  }

  delete scheme_;
}


real Block::get_t_begin() {

  if (parameter_.restart) {
    return restart::load_physical_time(parameter_.output_directory, block_id_.my_coords);
  } else {
    return parameter_.t_begin;
  }
}


real Block::get_t_end() { return parameter_.t_end; }


void Block::loop_cycle(real &t) {

  auto time_begin = std::chrono::steady_clock::now();


  dt_ = scheme_->get_dt();

  if (parameter_.criterion != nullptr && 
      fmod(t - parameter_.t_begin, dt_update_scheme_) < dt_) {
    update_scheme();
  }

  scheme_->step();

  t += dt_;


  // output
  parameter_.freeze_dt = freeze_dt_default_;
  int additional_outputs_before = parameter_.nconsecutive_outputs/2;

  if (fmod(t - parameter_.t_begin + additional_outputs_before*dt_, dt_output_vtk_) <
      parameter_.nconsecutive_outputs*dt_) { // vtk output
    scheme_->output(t, output_number_, true);
    ++output_number_;
    if (parameter_.nconsecutive_outputs > 1) {
      parameter_.freeze_dt = true;
    }
  } else if (fmod(t - parameter_.t_begin, dt_output_csv_) <
             parameter_.nconsecutive_outputs*dt_) { // no vtk output
    scheme_->output(t, output_number_, false);
  }


  auto time_end = std::chrono::steady_clock::now();

  if (nsteps_total >= 50) {
    wall_time_total_s += std::chrono::duration_cast<std::chrono::milliseconds>(time_end-time_begin).count()/1000.;
  }
  nsteps_total += 1;

  if (block_id_.my_id == 0) {
    std::cout << std::fixed << std::setprecision(5) << "\r\e[At = " << t
              << std::scientific << std::setprecision(3) << " / dt = " << dt_
              << std::fixed << std::setprecision(1)
              << " / progress: " << (t - parameter_.t_begin) * 100. / (parameter_.t_end - parameter_.t_begin) << " %";

    //std::cout << " / step duration/ms: "
    //          << std::chrono::duration_cast<std::chrono::milliseconds>(time_end-time_begin).count();

    //std::cout << std::scientific << std::setprecision(7) << " / avg step duration/s: "
    //          << wall_time_total_s/(nsteps_total-50);
    
    if (parameter_.electron_subcycling) {
      int electron_substeps = std::round(scheme_->get_dt(Parameter::kIon)/scheme_->get_dt());
      std::cout << " / electron substeps: " << electron_substeps;
    }

    std::cout << " / latest output: " << output_number_ - 1 << "          " << std::endl
              << std::flush;
  }
  
  // comment in for scaling study and adjust number of maximum steps as desired
  // //if (nsteps_total > 150) std::exit(EXIT_SUCCESS);
}


void Block::update_scheme() {

  block_id_.my_scheme_id = scheme_->evaluate_criterion();

  int requested_scheme_id[Parameter::kNDirections],
      requested_scheme_hierarchy_index[Parameter::kNDirections],
      my_scheme_hierarchy_index;

  boundary_.send_my_id(block_id_.my_scheme_id, block_id_.neighbour_scheme_ids); // update neighbour_scheme_ids

  for (int i = 0; i < parameter_.nschemes-2; ++i) { // go through scheme hierarchy from top to bottom

    if (block_id_.my_scheme_id == parameter_.scheme_hierarchy[i]) {

      for (int d = 0; d < parameter_.kNDirections; ++d) {

        // skip compatibility check at potential non-periodic boundary
        if ((block_id_.my_coords[0] == 0 &&                     d == Parameter::kXL && parameter_.is_periodic[0] == false) ||
            (block_id_.my_coords[0] == parameter_.nproc[0]-1 && d == Parameter::kXU && parameter_.is_periodic[0] == false) ||
            (block_id_.my_coords[1] == 0 &&                     d == Parameter::kYL && parameter_.is_periodic[1] == false) ||
            (block_id_.my_coords[1] == parameter_.nproc[1]-1 && d == Parameter::kYU && parameter_.is_periodic[1] == false) ||
            (block_id_.my_coords[2] == 0 &&                     d == Parameter::kZL && parameter_.is_periodic[2] == false) ||
            (block_id_.my_coords[2] == parameter_.nproc[2]-1 && d == Parameter::kZU && parameter_.is_periodic[2] == false)) {
          continue;
        }

        if (block_id_.neighbour_scheme_ids[d] != parameter_.scheme_hierarchy[std::max(0,i-1)] &&
            block_id_.neighbour_scheme_ids[d] != parameter_.scheme_hierarchy[i] &&
            block_id_.neighbour_scheme_ids[d] != parameter_.scheme_hierarchy[i+1]) { // check if all neighbour schemes are compatible
          block_id_.neighbour_scheme_ids[d] = parameter_.scheme_hierarchy[i+1];      // if not, assign a compatible scheme id
        }
      }
    }

    boundary_.send_neighbour_ids(block_id_.neighbour_scheme_ids, requested_scheme_id); // send scheme requirements to neighbours
    // from neighbour requirements select the scheme that is highest in the hierarchy
    for (int d = 0; d < parameter_.kNDirections; ++d) {
      requested_scheme_hierarchy_index[d] = -1;
      for (int s = 0; s < parameter_.nschemes; ++s) {
        if (parameter_.scheme_hierarchy[s] == requested_scheme_id[d]) {
          requested_scheme_hierarchy_index[d] = s;
          break;
        }
      }
    }
    my_scheme_hierarchy_index = *(std::min_element(requested_scheme_hierarchy_index, requested_scheme_hierarchy_index + Parameter::kNDirections));
    block_id_.my_scheme_id = parameter_.scheme_hierarchy[my_scheme_hierarchy_index];
    boundary_.send_my_id(block_id_.my_scheme_id, block_id_.neighbour_scheme_ids); // update neighbour_scheme_ids
  }

  int scheme_changed = 0; // int instead of boolean for easier mpi exchange
  if (block_id_.my_scheme_id != scheme_->get_scheme_id()) {

    scheme::Scheme *scheme_tmp = criterion::allocate_scheme(block_id_.my_scheme_id,
                                 block_id_, boundary_, parameter_);

    criterion::convert_scheme(scheme_, scheme_tmp, parameter_);

    delete scheme_;
    scheme_ = scheme_tmp;

    scheme_changed = 1;
  }

  scheme_changed = boundary_.mpi_max(scheme_changed);

  if (parameter_.nslots_per_node > 0 && // if is set
      scheme_changed) {
    load_balancing();
  }
}


void Block::load_balancing() {

  int* gathered_scheme_ids = nullptr,
     * gathered_mpi_coords = nullptr,
     * gathered_neighbour_ranks = nullptr,
     * send_data_to = nullptr,
     * receive_data_from = nullptr,
     * new_scheme_ids = nullptr,
     * new_mpi_coords = nullptr,
     * new_neighbour_ranks = nullptr;
  const int nprocs = parameter_.nproc[0]*parameter_.nproc[1]*parameter_.nproc[2],
            nnodes = nprocs/(max_local_rank_+1);

  if (block_id_.my_id == 0) {
    gathered_scheme_ids = new int[nprocs]();
    gathered_mpi_coords = new int[3*nprocs]();
    gathered_neighbour_ranks = new int[parameter_.kNDirections*nprocs]();
    send_data_to = new int[nprocs]();
    receive_data_from = new int[nprocs]();
    new_scheme_ids = new int[nprocs]();
    new_mpi_coords = new int[3*nprocs]();
    new_neighbour_ranks = new int[parameter_.kNDirections*nprocs]();
  }

  boundary_.gather(gathered_scheme_ids, block_id_.my_scheme_id);
  boundary_.gather(gathered_mpi_coords, block_id_.my_coords, 3);
  boundary_.gather(gathered_neighbour_ranks, block_id_.neighbour_ids, 6);

  if (block_id_.my_id == 0) {
    int current_node = 0,
        processes_on_node[nnodes];

    for (int n = 0; n < nnodes; ++n) {
      processes_on_node[n] = 0;
    }

    for (int s = 0; s < parameter_.nschemes; ++s) { // go through scheme hierarchy from top to bottom

      for (int i = 0; i < nprocs; ++i) {

        if (gathered_scheme_ids[i] != parameter_.scheme_hierarchy[s]) {
          continue;
        }

        send_data_to[i] = current_node*(max_local_rank_+1) + processes_on_node[current_node];
        receive_data_from[send_data_to[i]] = i;
        ++processes_on_node[current_node];
        if (processes_on_node[current_node]%parameter_.nslots_per_node == 0) {
          ++current_node;
        }

        if (current_node == nnodes) {
          current_node = 0;
        }
      }
    }

    for (int i = 0; i < nprocs; ++i) {

      new_scheme_ids[i] = gathered_scheme_ids[receive_data_from[i]];

      for (int d = 0; d < 3; ++d) {
        new_mpi_coords[i*3 + d] = gathered_mpi_coords[receive_data_from[i]*3 + d];
      }

      for (int d = 0; d < parameter_.kNDirections; ++d) {
        new_neighbour_ranks[i*parameter_.kNDirections + d] =
            send_data_to[gathered_neighbour_ranks[receive_data_from[i]*parameter_.kNDirections + d]];
      }
    }
  }

  int my_send_data_to = -1,
      my_receive_data_from = -1,
      my_new_scheme_id = -1;

  boundary_.scatter(send_data_to, &my_send_data_to);
  boundary_.scatter(receive_data_from, &my_receive_data_from);
  boundary_.scatter(new_scheme_ids, &my_new_scheme_id);

  scheme::Scheme* scheme_tmp = criterion::allocate_scheme(my_new_scheme_id,
                                  block_id_, boundary_, parameter_);

  int out_buffer_size = scheme_->get_data_size(),
      in_buffer_size = scheme_tmp->get_data_size();

  if (out_buffer_size <= 0 || in_buffer_size <= 0) {
      std::cerr<<std::endl<<"Error: Zero buffer size in load balancing exchange. "<<
      "Check whether all used schemes provide an accurate implementation of get_data_size()."<<std::endl<<std::endl;
      std::exit(EXIT_FAILURE);
  }

  real* out_buffer = new real[out_buffer_size](),
      * in_buffer = new real[in_buffer_size]();

  scheme_->send_model_data(out_buffer, my_send_data_to);
  scheme_tmp->receive_model_data(in_buffer, my_receive_data_from);

  scheme_tmp->set_dt(scheme_->get_dt());

  if (parameter_.electron_subcycling) {
    scheme_tmp->set_substeps(scheme_->get_substeps());
    scheme_tmp->set_substep_counter(scheme_->get_substep_counter());
  }

  delete scheme_;
  scheme_ = scheme_tmp;

  boundary_.scatter(new_mpi_coords, block_id_.my_coords, 3);
  boundary_.scatter(new_neighbour_ranks, block_id_.neighbour_ids, 6);

  block_id_.my_scheme_id = my_new_scheme_id;
  boundary_.send_my_id(block_id_.my_scheme_id, block_id_.neighbour_scheme_ids);

  for (int i = 0; i < 3; ++i) {
    parameter_.xb_loc[i] = parameter_.xb[i] + block_id_.my_coords[i] * parameter_.res_x[i] * parameter_.dx[i];
  }

  // clean up
  delete[] in_buffer;
  delete[] out_buffer;

  if (block_id_.my_id == 0) {
    delete[] new_neighbour_ranks;
    delete[] new_mpi_coords;
    delete[] new_scheme_ids;
    delete[] receive_data_from;
    delete[] send_data_to;
    delete[] gathered_neighbour_ranks;
    delete[] gathered_mpi_coords;
    delete[] gathered_scheme_ids;
  }
}


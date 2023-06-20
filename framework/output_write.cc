/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "output_write.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "framework/definitions.h"
#include "framework/parameter.h"

namespace output {

  void write_vti_file(const real* fields[], const int nfields, const std::string field_names[],
            const int output_number, const std::string output_directory,
            const real xb_loc[], const real dx[], const int coords[], const int res_x[], const int bd[])
  {
    int index_lower[3] = {coords[0]*res_x[0],
                          coords[1]*res_x[1],
                          coords[2]*res_x[2]};
    int index_upper[3] = {index_lower[0]+res_x[0],
                          index_lower[1]+res_x[1],
                          index_lower[2]+res_x[2]};

    int ncells = res_x[0]*res_x[1]*res_x[2];
    int offset = 0;
    int binary_size = ncells*sizeof(real) + sizeof(int);

    std::string data_type = "Float32";
    if (sizeof(real) == 8)
      data_type = "Float64";

    std::string filepath = output_directory + "/vti/output" + std::to_string(output_number) +
        "_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) +
        "_" + std::to_string(coords[2]) + ".vti";

    if (coords[0] < 0) { //ugly fix due to blinkers in write routine
      index_lower[0] = index_lower[1] = index_lower[2] = 0;
      index_upper[0] = res_x[0]; index_upper[1] = res_x[1]; index_upper[2] = res_x[2];
      filepath = output_directory + "/output" + std::to_string(output_number) +
          "_" + std::to_string(-coords[0]) + "_" + std::to_string(-coords[1]) +
          "_" + std::to_string(-coords[2]) + ".vti";
    }

    // header
    std::ofstream os(filepath.c_str());

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<ImageData WholeExtent=\"" << index_lower[0] << " " << index_upper[0] << " " << index_lower[1] << " "
       << index_upper[1] << " " << index_lower[2] << " " << index_upper[2] << "\" Origin=\"" << xb_loc[0] << " "
       << xb_loc[1] << " " << xb_loc[2] << "\" Spacing=\"" << dx[0] << " " << dx[1] << " " << dx[2] << "\">\n";
    os << "<Piece Extent=\"" << index_lower[0] << " " << index_upper[0] << " " << index_lower[1] << " "
       << index_upper[1] << " " << index_lower[2] << " " << index_upper[2] << "\">\n";

    os << "<CellData Scalars=\""<< "output" << "\">\n";

    for (int ifield = 0; ifield < nfields; ++ifield)
    {
      os << "<DataArray type=\"" << data_type << "\" Name=\"" << field_names[ifield] <<
          "\" format=\"appended\" offset=\"" << offset << "\"/>\n";
      offset += binary_size;
    }

    os << "</CellData>\n";
    os << "</Piece>\n";
    os << "</ImageData>\n";
    os << "<AppendedData encoding=\"raw\">\n";
    os << "_" ;

    os.close();

    // binary data
    int N = ncells*sizeof(real);
    os.open(filepath.c_str(), std::ios::binary | std::ios::app);

    for (int ifield = 0; ifield < nfields; ++ifield)
    {
      os.write(reinterpret_cast<const char*>(&N), sizeof(int));
      int i = 0;
      for (int z = -bd[2]; z < res_x[2]+bd[2]; ++z)
      {
        for (int y = -bd[1]; y < res_x[1]+bd[1]; ++y)
        {
          for (int x = -bd[0]; x < res_x[0]+bd[0]; ++x)
          {
            // boundary cells are skipped
            if (x >= 0 && x < res_x[0] && y >= 0 && y < res_x[1] && z >= 0 && z < res_x[2])
            {
              os.write(reinterpret_cast<const char*>(&fields[ifield][i]), sizeof(real));
            }
            ++i;
          }
        }
      }
    }

    os.close();

    // footer
    os.open(filepath.c_str(), std::ios::app);

    os << "\n</AppendedData>\n";
    os << "</VTKFile>\n";

    os.close();
  }


  void write_pvti_file(const int nfields, const std::string field_names[], const real physical_time,
            const int output_number, const std::string output_directory,
            const real xb[], const real dx[], const int res_x[], const int res_x_total[])
  {
    std::string filepath = output_directory + "/output" + std::to_string(output_number) + ".pvti";
    std::ofstream os(filepath.c_str());

    std::string data_type = "Float32";
    if (sizeof(real) == 8)
      data_type = "Float64";

    os << "<?xml version=\"1.0\" ?>\n";
    os << "<!-- t = " << std::fixed << std::setprecision(15) << physical_time <<" -->\n";
    os << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PImageData WholeExtent=\"" << 0 << " " << res_x_total[0] << " " << 0 << " "
    << res_x_total[1] << " " << 0 << " " << res_x_total[2] << "\" GhostLevel=\"0\" Origin=\"" <<
    xb[0] << " " <<  xb[1] << " " << xb[2] << "\" Spacing=\"" << dx[0] << " "<< dx[1] << " "<< dx[2] << "\">\n";
    os << "<PCellData Scalars=\"" << "output" << "\">\n";

    for (int ifield = 0; ifield < nfields; ++ifield)
    {
      os << "<PDataArray type=\"" << data_type << "\" Name=\"" << field_names[ifield] << "\"/>\n";
    }

    os << "</PCellData>\n";

    for (int pos_x = res_x[0]; pos_x < res_x_total[0]+1; pos_x += res_x[0])
    {
      for (int pos_y = res_x[1]; pos_y < res_x_total[1]+1; pos_y += res_x[1])
      {
        for (int pos_z = res_x[2]; pos_z < res_x_total[2]+1; pos_z += res_x[2])
        {
          os << "<Piece Extent=\"" << pos_x - res_x[0] << " " << pos_x << " " << pos_y - res_x[1]
            << " " << pos_y << " " << pos_z - res_x[2] << " " << pos_z <<"\" Source=\""
            << "vti/" << "output" << output_number << "_"  << pos_x/res_x[0] - 1 << "_"
            << pos_y/res_x[1] - 1 << "_" << pos_z/res_x[2] - 1 << ".vti\"/>\n";

        }
      }
    }

    os << "</PImageData>\n";
    os << "</VTKFile>\n";

    os.close();
  }


  void write_text_file(const real fields[], const int nfields, const std::string field_names[],
            const int output_number, const std::string output_directory, const int coords[])
  {
    std::string filepath = output_directory + "/csv/output" +
      "_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) +
      "_" + std::to_string(coords[2]) + ".csv";

    std::ofstream os;

    if (output_number == 0) {
      // header
      os.open(filepath.c_str());
      for (int ifield = 0; ifield < nfields; ++ifield) {
        os << "# (" << ifield << ") " << field_names[ifield] << "\n";
      }
      os.close();
    } else {
      os.close();
      os.open(filepath.c_str(), std::ios::app);
    }
    // data
    for (int ifield = 0; ifield < nfields-1; ++ifield) {
      os << fields[ifield] << " ";
    }
    os << fields[nfields-1] << "\n";

    os.close();
  }


  void write_restart_file(const int scheme_id, const real physical_time,
                          const real scalars[], const int nscalars,
                          const real* fields[], const int field_sizes[], const int nfields,
                          const std::string output_directory, const int coords[],
                          const int output_number) /* output_number = -1 */
  {
    std::string filepath = output_directory + "/restart/restart_" +
        std::to_string(coords[0]) + "_" + std::to_string(coords[1]) + "_" +
        std::to_string(coords[2]) + ".bin";

    if( output_number > 1 ) {
      std::string newfilepath = output_directory + "/phase_space/restart_" +
          std::to_string(coords[0]) + "_" + std::to_string(coords[1]) + "_" +
          std::to_string(coords[2]) + ".bin";
      std::string copycall = "mv " + filepath + ".old " + newfilepath + "." +
          std::to_string(output_number-2) + " 2> /dev/null";
      system(copycall.c_str());
    }

    // TODO this is suboptimal
    system(("mv " + filepath + " " + filepath + ".old 2> /dev/null").c_str());

    std::ofstream os(filepath.c_str(), std::ios::binary | std::ios::out);

    os.write(reinterpret_cast<const char*>(&scheme_id), sizeof(int));
    os.write(reinterpret_cast<const char*>(&physical_time), sizeof(real));
    os.write(reinterpret_cast<const char*>(&nscalars), sizeof(int));
    os.write(reinterpret_cast<const char*>(&scalars[0]), nscalars*sizeof(real));

    for (int i = 0; i < nfields; ++i) {
      os.write(reinterpret_cast<const char*>(&field_sizes[i]), sizeof(int));
      os.write(reinterpret_cast<const char*>(&fields[i][0]), field_sizes[i]*sizeof(real));
    }

    os.close();
  }

} // namespace output


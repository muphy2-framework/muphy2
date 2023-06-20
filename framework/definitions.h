/*
 * This file is part of the muphyII multiphysics plasma simulation project.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#ifdef SINGLE_PRECISION
    typedef float real;
    #define MPI_REAL_CXX MPI_FLOAT

#else
    typedef double real;
    #define MPI_REAL_CXX MPI_DOUBLE

#endif // SINGLE_PRECISION

#endif // DEFINITIONS_H_

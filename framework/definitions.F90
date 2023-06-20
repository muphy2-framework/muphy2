
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef DEFINITIONS_F90_
#define DEFINITIONS_F90_

#ifdef SINGLE_PRECISION
#define SELECTED_PRECISION REAL32
#define MPI_REALTYPE MPI_REAL4

#else
#define SELECTED_PRECISION REAL64
#define MPI_REALTYPE MPI_REAL8

#endif
! SINGLE_PRECISION
#endif
! DEFINITIONS_F90_

module definitions
    use, intrinsic :: iso_fortran_env
    implicit none

    integer, parameter :: PRC = SELECTED_PRECISION
    real(kind=PRC), parameter :: M_PI = 3.14159265358979323846264338327950288
end module 

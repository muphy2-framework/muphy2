
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

subroutine output_extract_distribution_function_slice(f_slice,f,dimX,dimV,max_dimV,BD,direction,slice_coordinates)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), max_dimV, BD(3), direction, slice_coordinates(4)
    real(kind=PRC), intent(out), dimension(-BD(direction):dimX(direction)+BD(direction),max_dimV) :: f_slice
    real(kind=PRC), intent(in ), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),dimV(1),dimV(2),dimV(3)) :: f

    integer :: offset

    offset = (max_dimV-dimV(direction))/2

    if (direction == 1) then
        f_slice(:,1+offset:offset+dimV(direction)) = f(:,slice_coordinates(1),slice_coordinates(2),&
                    :,slice_coordinates(3)+1,slice_coordinates(4)+1) ! +1 in velocity coordinates because array starts at 1
    else if (direction == 2) then
        f_slice(:,1+offset:offset+dimV(direction)) = f(slice_coordinates(1),:,slice_coordinates(2),&
                    slice_coordinates(3)+1,:,slice_coordinates(4)+1)
    else if (direction == 3) then
        f_slice(:,1+offset:offset+dimV(direction)) = f(slice_coordinates(1),slice_coordinates(2),:,&
                    slice_coordinates(3)+1,slice_coordinates(4)+1,:)
    end if

end subroutine output_extract_distribution_function_slice


subroutine output_calc_statistics_fluid(n,u,T,m,mass,kinetic_energy,thermal_energy,&
                momentum_x,momentum_y,momentum_z,dimX,BD,dx,dimensionality_x,dimensionality_v)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_x, dimensionality_v
    real(kind=PRC), intent(in) :: dx(3), m
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: n, T
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u
    real(kind=PRC), intent(out) :: mass, kinetic_energy, thermal_energy, momentum_x, momentum_y, momentum_z

    real(kind=PRC) :: d3

    d3 = dx(1)
    if (dimensionality_x > 1) d3 = d3*dx(2)
    if (dimensionality_x > 2) d3 = d3*dx(3)

    mass = m*sum(n(0:dimX(1),0:dimX(2),0:dimX(3)))*d3
    momentum_x = m*sum(u(0:dimX(1),0:dimX(2),0:dimX(3),1)*n(0:dimX(1),0:dimX(2),0:dimX(3)))*d3
    momentum_y = m*sum(u(0:dimX(1),0:dimX(2),0:dimX(3),2)*n(0:dimX(1),0:dimX(2),0:dimX(3)))*d3
    momentum_z = m*sum(u(0:dimX(1),0:dimX(2),0:dimX(3),3)*n(0:dimX(1),0:dimX(2),0:dimX(3)))*d3
    kinetic_energy = .5_PRC*m*sum(&
                     (u(0:dimX(1),0:dimX(2),0:dimX(3),1)**2 + &
                      u(0:dimX(1),0:dimX(2),0:dimX(3),2)**2 + &
                      u(0:dimX(1),0:dimX(2),0:dimX(3),3)**2) * &
                      n(0:dimX(1),0:dimX(2),0:dimX(3)) &
                      ) * d3
    thermal_energy = .5_PRC*dimensionality_v*sum(T(0:dimX(1),0:dimX(2),0:dimX(3))*&
                         n(0:dimX(1),0:dimX(2),0:dimX(3)))*d3

end subroutine output_calc_statistics_fluid


subroutine output_calc_statistics_electromagnetic(E,B,eps0,mu0,electric_energy,magnetic_energy,&
                          em_momentum_x,em_momentum_y,em_momentum_z,dimX,BD,dx,dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_x
    real(kind=PRC), intent(in) :: dx(3), eps0, mu0
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E, B
    real(kind=PRC), intent(out) :: electric_energy, magnetic_energy, em_momentum_x, em_momentum_y, em_momentum_z

    real(kind=PRC) :: d3

    d3 = dx(1)
    if (dimensionality_x > 1) d3 = d3*dx(2)
    if (dimensionality_x > 2) d3 = d3*dx(3)

    electric_energy = .5_PRC*eps0*sum(E(0:dimX(1),0:dimX(2),0:dimX(3),1)**2 + &
                                      E(0:dimX(1),0:dimX(2),0:dimX(3),2)**2 + &
                                      E(0:dimX(1),0:dimX(2),0:dimX(3),3)**2) * d3
    magnetic_energy = .5_PRC/mu0 *sum(B(0:dimX(1),0:dimX(2),0:dimX(3),1)**2 + &
                                      B(0:dimX(1),0:dimX(2),0:dimX(3),2)**2 + &
                                      B(0:dimX(1),0:dimX(2),0:dimX(3),3)**2) * d3
    em_momentum_x = eps0*sum(E(0:dimX(1),0:dimX(2),0:dimX(3),2)*B(0:dimX(1),0:dimX(2),0:dimX(3),3) - &
                             E(0:dimX(1),0:dimX(2),0:dimX(3),3)*B(0:dimX(1),0:dimX(2),0:dimX(3),2))*d3
    em_momentum_y = eps0*sum(E(0:dimX(1),0:dimX(2),0:dimX(3),3)*B(0:dimX(1),0:dimX(2),0:dimX(3),1) - &
                             E(0:dimX(1),0:dimX(2),0:dimX(3),1)*B(0:dimX(1),0:dimX(2),0:dimX(3),3))*d3
    em_momentum_z = eps0*sum(E(0:dimX(1),0:dimX(2),0:dimX(3),1)*B(0:dimX(1),0:dimX(2),0:dimX(3),2) - &
                             E(0:dimX(1),0:dimX(2),0:dimX(3),2)*B(0:dimX(1),0:dimX(2),0:dimX(3),1))*d3

end subroutine output_calc_statistics_electromagnetic

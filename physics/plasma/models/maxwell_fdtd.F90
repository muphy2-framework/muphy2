
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.

! Finite Difference Time Domain (FDTD) method
!
! _______ Stepping (Leapfrog) and Order of the Scheme _______
! The scheme is second order accurate due to the leapfrog-like
! stepping. In order to set up the necessary displacement of B in
! time by .5 dt, there is an Euler init step. TODO: It makes
! no difference whether the init step is done or not.
! 
! The general stepping is:
! maxwell init -> loop( .5 maxwell -> plasma -> .5 maxwell ), so that
! the fdtd receives the correct j on average.
! In this implementation, the first half maxwell step is included in
! the init step and two half maxwell steps in the loop are combined into
! one full step.
!
! _______ Subcycling _______
! The dt_maxwell needs to be constant throughout the simulation so actually
! the plasma step is a supercycle over multiple maxwell steps, i.e. dt_plasma
! is always a multiple of dt_maxwell. In fluid simulations the dt_plasma can change 
! which is implemented by a change of the number of maxwell subcycles.
! dt_maxwell must be chosen small enough so that dt_plasma does not become
! smaller. Therefore a safety factor is introduced for setting the initial
! dt_maxwell based on dt_plasma. Often it is more efficient to choose a small
! dt_maxwell so that multiples of it are close to the actual dt_plasma.
!
! Because of the leapfrog stepping, changing between odd and even numbers
! of maxwell subcycles requires a change from B at halftime to E at halftime
! (or the other way around). This is possible but adds a lot of code and an
! array for E at halftime. Thus subcycles are restricted to odd values.

subroutine step_E_maxwell_fdtd(E, B, j, c, mu0, dimX, BD, dx, dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(1:3), BD(3)
    real(kind=PRC), intent(in) :: dx(3), dt, c, mu0
    real(kind=PRC), intent(in),    dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: j, B
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E

    integer :: x, y, z
    real(kind=PRC) :: dtcc
    dtcc = dt*c*c

    !$acc data present(E, B, j)

    !$acc parallel
    !$acc loop gang
    do z = 0,dimX(3)
        !$acc loop vector collapse(2)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                E(x,y,z,1) = E(x,y,z,1) + dtcc*( &
                        (B(x,y+1,z,3) - B(x,y,z,3))/dx(2) - &
                        (B(x,y,z+1,2) - B(x,y,z,2))/dx(3) - &
                        j(x,y,z,1)*mu0 )
                E(x,y,z,2) = E(x,y,z,2) + dtcc*( &
                        (B(x,y,z+1,1) - B(x,y,z,1))/dx(3) - &
                        (B(x+1,y,z,3) - B(x,y,z,3))/dx(1) - &
                        j(x,y,z,2)*mu0 )
                E(x,y,z,3) = E(x,y,z,3) + dtcc*( &
                        (B(x+1,y,z,2) - B(x,y,z,2))/dx(1) - &
                        (B(x,y+1,z,1) - B(x,y,z,1))/dx(2) - &
                        j(x,y,z,3)*mu0 )
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

end subroutine step_E_maxwell_fdtd


subroutine step_B_ipol_maxwell_fdtd(E, B, B_ipol, dimX, BD, dx, dt) ! with interpolation to integer times
    use definitions
    implicit none

    integer, intent(in) :: dimX(1:3), BD(3)
    real(kind=PRC), intent(in) :: dx(1:3), dt
    real(kind=PRC), intent(in),    dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: E
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: B
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: B_ipol

    integer :: x, y, z

    !$acc data present(E, B, B_ipol)

    !$acc parallel
    !$acc loop gang
    do z = 0,dimX(3)
        !$acc loop vector collapse(2)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                ! B_ipol = .5*B_old
                B_ipol(x,y,z,1) = 0.5_PRC * B(x,y,z,1)
                B_ipol(x,y,z,2) = 0.5_PRC * B(x,y,z,2)
                B_ipol(x,y,z,3) = 0.5_PRC * B(x,y,z,3)
            end do
        end do
    end do
    !$acc end parallel

    call step_B_maxwell_fdtd(E, B, dimX, BD, dx, dt)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                ! B_ipol += .5*B_new
                B_ipol(x,y,z,1) = B_ipol(x,y,z,1) + 0.5_PRC*B(x,y,z,1)
                B_ipol(x,y,z,2) = B_ipol(x,y,z,2) + 0.5_PRC*B(x,y,z,2)
                B_ipol(x,y,z,3) = B_ipol(x,y,z,3) + 0.5_PRC*B(x,y,z,3)
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

end subroutine step_B_ipol_maxwell_fdtd


subroutine step_B_maxwell_fdtd(E, B, dimX, BD, dx, dt)
    use definitions
    implicit none

    integer, intent(in) :: dimX(1:3), BD(3)
    real(kind=PRC), intent(in) :: dx(1:3), dt
    real(kind=PRC), intent(in)    :: E(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)
    real(kind=PRC), intent(inout) :: B(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3)

    integer :: x, y, z

    !$acc data present(E, B)

    !$acc parallel
    !$acc loop gang
    do z = 0,dimX(3)
        !$acc loop vector collapse(2)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                B(x,y,z,1) = B(x,y,z,1) - dt * (&
                        (E(x,y,z,3) - E(x,y-1,z,3))/dx(2) - &
                        (E(x,y,z,2) - E(x,y,z-1,2))/dx(3) )
                B(x,y,z,2) = B(x,y,z,2) - dt * ( &
                        (E(x,y,z,1) - E(x,y,z-1,1))/dx(3) - &
                        (E(x,y,z,3) - E(x-1,y,z,3))/dx(1) )
                B(x,y,z,3) = B(x,y,z,3) - dt * ( &
                        (E(x,y,z,2) - E(x-1,y,z,2))/dx(1) - &
                        (E(x,y,z,1) - E(x,y-1,z,1))/dx(2) )
            end do
        end do
    end do
    !$acc end parallel

    !$acc end data

end subroutine step_B_maxwell_fdtd


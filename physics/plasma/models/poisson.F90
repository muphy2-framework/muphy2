
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


! Solve for given rho:
! laplace(phi) = -rho / eps0
! E = -grad(phi)


subroutine poisson_calc_E_from_phi(E, phi, dimX, BD, dx) ! E = -grad(phi)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in) :: dx(3)
    real(kind=PRC), intent(in   ), dimension(-BD(1):dimX(1)+BD(1), -BD(2):dimX(2)+BD(2), -BD(3):dimX(3)+BD(3)) :: phi
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1), -BD(2):dimX(2)+BD(2), -BD(3):dimX(3)+BD(3), 3) :: E

    integer :: x, y, z

    !$acc data present(E, phi)

    !$acc parallel
    !$acc loop gang
    do z = 0,dimX(3)
        !$acc loop vector collapse(2)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                E(x,y,z,1) = (-phi(x+1,y,z) + phi(x-1,y,z)) / (2._PRC*dx(1))
                E(x,y,z,2) = (-phi(x,y+1,z) + phi(x,y-1,z)) / (2._PRC*dx(2))
                E(x,y,z,3) = (-phi(x,y,z+1) + phi(x,y,z-1)) / (2._PRC*dx(3))
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine poisson_calc_E_from_phi


subroutine poisson_gauss_seidel(phi, rho, eps0, dimX, BD, dx, dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_x
    real(kind=PRC), intent(in) :: dx(3), eps0
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: phi
    real(kind=PRC), intent(in   ), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: rho

    integer :: x, y, z, i
    real(kind=PRC) :: invdx2(3), suminv

    ! this does not parallelize on GPUs

    invdx2 = 1._PRC/(dx*dx)
    suminv = 0._PRC
    do i = 1,dimensionality_x
        suminv = suminv + 2._PRC*invdx2(i)
    end do

    !$acc data present(phi, rho) copyin(invdx2)

    if (dimensionality_x == 1) then
        !$acc kernels
        do z=-BD(3)+1,dimX(3)+BD(3)-1
            do y=-BD(2)+1,dimX(2)+BD(2)-1
                do x=-BD(1)+1,dimX(1)+BD(1)-1
                    phi(x,y,z) = ((phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                                   rho(x,y,z)/eps0) / suminv
                enddo
            enddo
        enddo
        !$acc end kernels
    else if (dimensionality_x == 2) then
        !$acc kernels
        do z=-BD(3)+1,dimX(3)+BD(3)-1
            do y=-BD(2)+1,dimX(2)+BD(2)-1
                do x=-BD(1)+1,dimX(1)+BD(1)-1
                    phi(x,y,z) = ((phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                                  (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                                   rho(x,y,z)/eps0) / suminv
                enddo
            enddo
        enddo
        !$acc end kernels
    else ! dimensionality_x == 3
        !$acc kernels
        do z=-BD(3)+1,dimX(3)+BD(3)-1
            do y=-BD(2)+1,dimX(2)+BD(2)-1
                do x=-BD(1)+1,dimX(1)+BD(1)-1
                    phi(x,y,z) = ((phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                                  (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                                  (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                                   rho(x,y,z)/eps0) / suminv
                enddo
            enddo
        enddo
        !$acc end kernels
    end if

    !$acc end data

end subroutine poisson_gauss_seidel


subroutine poisson_sor(phi, rho, eps0, dimX, BD, dx, dimensionality_x)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3), dimensionality_x
    real(kind=PRC), intent(in) :: dx(3), eps0
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: phi
    real(kind=PRC), intent(in   ), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: rho

    integer :: x, y, z, i
    real(kind=PRC) :: omega, invdx2(3), suminv

    ! Successive Over-Relaxation (SOR)
    ! equivalent to red-black Gauss-Seidel for omega = 1

    omega = 2._PRC/(1._PRC+sin(M_PI/(dimX(1)+1)))
    suminv = 0._PRC
    invdx2 = 0._PRC
    do i = 1,dimensionality_x
        invdx2(i) = 1._PRC/(dx(i)*dx(i))
        suminv = suminv + 2._PRC*invdx2(i)
    end do

    !$acc data present(phi, rho) copyin(invdx2)

    !$acc parallel

    ! even cells
    !$acc loop gang
    do z=-BD(3)+2,dimX(3)+BD(3)-1,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+1,dimX(2)+BD(2)-2,2
            do x=-BD(1)+1,dimX(1)+BD(1)-2,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc loop gang
    do z=-BD(3)+1,dimX(3)+BD(3)-2,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+2,dimX(2)+BD(2)-1,2
            do x=-BD(1)+1,dimX(1)+BD(1)-2,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc loop gang
    do z=-BD(3)+1,dimX(3)+BD(3)-2,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+1,dimX(2)+BD(2)-2,2
            do x=-BD(1)+2,dimX(1)+BD(1)-1,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc loop gang
    do z=-BD(3)+2,dimX(3)+BD(3)-1,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+2,dimX(2)+BD(2)-1,2
            do x=-BD(1)+2,dimX(1)+BD(1)-1,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    ! odd cells
    !$acc loop gang
    do z=-BD(3)+1,dimX(3)+BD(3)-2,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+1,dimX(2)+BD(2)-2,2
            do x=-BD(1)+1,dimX(1)+BD(1)-2,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc loop gang
    do z=-BD(3)+1,dimX(3)+BD(3)-2,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+2,dimX(2)+BD(2)-1,2
            do x=-BD(1)+2,dimX(1)+BD(1)-1,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc loop gang
    do z=-BD(3)+2,dimX(3)+BD(3)-1,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+1,dimX(2)+BD(2)-2,2
            do x=-BD(1)+2,dimX(1)+BD(1)-1,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc loop gang
    do z=-BD(3)+2,dimX(3)+BD(3)-1,2
        !$acc loop vector collapse(2)
        do y=-BD(2)+2,dimX(2)+BD(2)-1,2
            do x=-BD(1)+1,dimX(1)+BD(1)-2,2
                phi(x,y,z) = phi(x,y,z) + omega*(( &
                    (phi(x+1,y,z) + phi(x-1,y,z))*invdx2(1) + &
                    (phi(x,y+1,z) + phi(x,y-1,z))*invdx2(2) + &
                    (phi(x,y,z+1) + phi(x,y,z-1))*invdx2(3) + &
                    rho(x,y,z)/eps0)/suminv - phi(x,y,z))
            enddo
        enddo
    enddo

    !$acc end parallel

    !$acc end data

end subroutine poisson_sor


function poisson_difference_max_value(phi_old, phi_new, dimX, BD) result(phi_diff_max)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: phi_old, phi_new
    real(kind=PRC) :: phi_diff_max

    !$acc data present(phi_old, phi_new) copyout(phi_diff_max)
    !$acc kernels
    phi_diff_max = maxval(abs(phi_new-phi_old))
    !$acc end kernels
    !$acc end data

end function poisson_difference_max_value


function poisson_sum_phi(phi, dimX, BD) result(sum_inner_cells)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3)) :: phi
    real(kind=PRC) :: sum_inner_cells
    integer :: x, y, z

    !sum_inner_cells = sum(phi(0:dimX(1),0:dimX(2),0:dimX(3)))

    sum_inner_cells = 0._PRC

    !$acc data present(phi) copy(sum_inner_cells)
    !$acc kernels
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            do x = 0,dimX(1)
                sum_inner_cells = sum_inner_cells + phi(x,y,z)
            enddo
        enddo
    enddo
    !$acc end kernels
    !$acc end data

end function poisson_sum_phi


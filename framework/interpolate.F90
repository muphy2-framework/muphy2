
! This file is part of the muphyII multiphysics plasma simulation project.
! 
! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at https://mozilla.org/MPL/2.0/.


subroutine interpolate_yee_face_to_centered(u_face,u_centered,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_face
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_centered

    integer :: x,y,z

    !$acc data present_or_copy(u_centered, u_face)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                u_centered(x,y,z,1) = 0.5_PRC*(u_face(x,y,z,1) + u_face(x+1,y,z,1))
                u_centered(x,y,z,2) = 0.5_PRC*(u_face(x,y,z,2) + u_face(x,y+1,z,2))
                u_centered(x,y,z,3) = 0.5_PRC*(u_face(x,y,z,3) + u_face(x,y,z+1,3))
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine interpolate_yee_face_to_centered


subroutine interpolate_yee_edge_to_centered(u_edge,u_centered,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_edge
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_centered

    integer :: x,y,z

    !$acc data present_or_copy(u_centered, u_edge)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                u_centered(x,y,z,1) = 0.25_PRC*(u_edge(x,y,z,1) + u_edge(x,y+1,z,1) + u_edge(x,y,z+1,1) + u_edge(x,y+1,z+1,1))
                u_centered(x,y,z,2) = 0.25_PRC*(u_edge(x,y,z,2) + u_edge(x+1,y,z,2) + u_edge(x,y,z+1,2) + u_edge(x+1,y,z+1,2))
                u_centered(x,y,z,3) = 0.25_PRC*(u_edge(x,y,z,3) + u_edge(x+1,y,z,3) + u_edge(x,y+1,z,3) + u_edge(x+1,y+1,z,3))
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine interpolate_yee_edge_to_centered


subroutine interpolate_yee_centered_to_face(u_centered,u_face,dimX,BD)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), BD(3)
    real(kind=PRC), intent(in), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_centered
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3),3) :: u_face

    integer :: x,y,z

    !$acc data present_or_copy(u_centered, u_face)

    !$acc parallel
    !$acc loop gang collapse(2)
    do z = 0,dimX(3)
        do y = 0,dimX(2)
            !$acc loop vector
            do x = 0,dimX(1)
                u_face(x,y,z,1) = 0.5_PRC*(u_centered(x,y,z,1) + u_centered(x-1,y,z,1))
                u_face(x,y,z,2) = 0.5_PRC*(u_centered(x,y,z,2) + u_centered(x,y-1,z,2))
                u_face(x,y,z,3) = 0.5_PRC*(u_centered(x,y,z,3) + u_centered(x,y,z-1,3))
            enddo
        enddo
    enddo
    !$acc end parallel

    !$acc end data

end subroutine interpolate_yee_centered_to_face


subroutine interpolate_linear_extrapolation_to_boundary_f(f,dimX,dimV,BD,direction)
    use definitions
    implicit none

    integer, intent(in) :: dimX(3), dimV(3), BD(3), direction
    real(kind=PRC), intent(inout), dimension(-BD(1):dimX(1)+BD(1),-BD(2):dimX(2)+BD(2),-BD(3):dimX(3)+BD(3), &
                                              dimV(1),dimV(2),dimV(3)) :: f
    integer :: i
    
    !$acc data present(f)

    select case (direction)
        case (0) ! x lower
            do i=1,BD(1)
                !$acc kernels
                f(-i,:,:,:,:,:) = f(0,:,:,:,:,:) + i*(f(0,:,:,:,:,:) - f(1,:,:,:,:,:))
                !$acc end kernels
            end do
        case (1) ! x upper
            do i=1,BD(1)
                !$acc kernels
                f(dimX(1)+i,:,:,:,:,:) = f(dimX(1),:,:,:,:,:) + i*(f(dimX(1),:,:,:,:,:) - f(dimX(1)-1,:,:,:,:,:))
                !$acc end kernels
            end do
        case (2) ! y lower
            do i=1,BD(2)
                !$acc kernels
                f(:,-i,:,:,:,:) = f(:,0,:,:,:,:) + i*(f(:,0,:,:,:,:) - f(:,1,:,:,:,:))
                !$acc end kernels
            end do
        case (3) ! y upper
            do i=1,BD(2)
                !$acc kernels
                f(:,dimX(2)+i,:,:,:,:) = f(:,dimX(2),:,:,:,:) + i*(f(:,dimX(2),:,:,:,:) - f(:,dimX(2)-1,:,:,:,:))
                !$acc end kernels
            end do
        case (4) ! z lower
            do i=1,BD(3)
                !$acc kernels
                f(:,:,-i,:,:,:) = f(:,:,0,:,:,:) + i*(f(:,:,0,:,:,:) - f(:,:,1,:,:,:))
                !$acc end kernels
            end do
        case (5) ! z upper
            do i=1,BD(3)
                !$acc kernels
                f(:,:,dimX(3)+i,:,:,:) = f(:,:,dimX(3),:,:,:) + i*(f(:,:,dimX(3),:,:,:) - f(:,:,dimX(3)-1,:,:,:))
                !$acc end kernels
            end do
    end select

    !$acc end data

end subroutine interpolate_linear_extrapolation_to_boundary_f
